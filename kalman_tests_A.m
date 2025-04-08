% 功能：在目标数量显著超过传感器数量时，根据信息增益（KL偏差）和风险来比较传感器管理者的性能。


% if a movie is to be made
make_movie = false;%true;

if make_movie
    aviobj = avifile('risk_test_multi_class.avi', 'fps', 5,'compression','none');
end
frames=0;

% update_order = [1:10]';

% calculate coordinates for a circle w/ given radius (used to draw sensor FOV)
radius = 5e2;
a=0:359;
cx=radius*cosd(a);
cy=radius*sind(a);

% calculate coordinates for a square w/ given radius (used to draw sensor
% FOV)
sx = [-radius radius radius -radius];
sy = [-radius -radius radius radius];


% load in target truth
% Px = load('scenario1_x');
% Py = load('scenario1_y');
% Vx = load('scenario1_vx');
% Vy = load('scenario1_vy');

Px = load('scenario2_x');
Py = load('scenario2_y');
Vx = load('scenario2_vx');
Vy = load('scenario2_vy');


% classification truth used in generating measurements
%class_truthc = [1 3 1 1 1 5 1 2 4 2];  % 3 toi
%class_truthc = [3 3 3 3 1 5 1 2 4 2];  % 6 toi
%class_truthc = [1 1 1 3 1 5 1 1 4 1];  % 1 toi
%class_truthc = [2 1 2 1 1 5 1 1 4 2];  % 3 toi 
class_truthc = [1 1 1 2 1 2 1 2 4 1];  % 3 toi
%class_truthc = [2 1 3 3 1 1 1 1 1 3];  % 1 toi & 3 mtoi
%class_truthc = [3 1 3 3 1 1 1 1 1 2];  % 1 toi & 3 mtoi
%class_truthc = [1 1 1 4 1 3 1 1 1 2];  % 1 toi & 1 mtoi & 1 ltoi
%class_truthc = [1 2 4 1 1 4 1 1 2 1];  % 2 toi & 2 ltoi
%class_truthc = [2 1 2 3 1 3 1 3 1 2];  % 3 toi & 3 mtoi
%class_truthc = [2 1 3 1 2 5 1 1 4 2];  % 4 toi

% how often to perform an update if periodically updating
time_step = 2;
periodic_update = false;
% periodic_update = true;

tau = 1;%.01;%.1;
kmax = 200;
mc_runs = 1; %100;  % monte carlo runs (set to 1 for video / display)
num_targs = length(class_truthc);

% keep track of root mean squared error 跟踪均方根误差
mse = zeros(num_targs, kmax * mc_runs);

% keep histogram of times each target is measured 保留每个目标测量时间的直方图
hgram = zeros(num_targs, kmax * mc_runs);

% keep track of ERR at each time step 跟踪每个时间步骤的ERR
err_collection = zeros(num_targs, kmax / time_step);

% the confusion matrix for each sensor.
cm1 = [.80 .05 .05 .05 .05];
cm2 = [.05 .80 .05 .05 .05];
cm3 = [.05 .05 .80 .05 .05];
cm4 = [.05 .05 .05 .80 .05];
cm5 = [.05 .05 .05 .05 .80];
cm = [cm1;cm2;cm3;cm4;cm5];

% cm1 = [.95 .05/4 .05/4 .05/4 .05/4];
% cm2 = [.05/4 .95 .05/4 .05/4 .05/4];
% cm3 = [.05/4 .05/4 .95 .05/4 .05/4];
% cm4 = [.05/4 .05/4 .05/4 .95 .05/4];
% cm5 = [.05/4 .05/4 .05/4 .05/4 .95];
% cm = [cm1;cm2;cm3;cm4;cm5];

% initial states for all sensors
s = [0 0 0 0]';

% linear estimate to measurement state transform and measurement error
% covariance
%H = diag([1 1 1 1]);
%R = diag([30 1 30 1]);

%Variance of noise process that models velocity once integrated
Sv = 5;%1;%5;

%State transition matrix for the state vector = [Px Vx Py Vy]'
F = [1 time_step; 0 1];
F = blkdiag(F,F);

%Process noise covariance matrix for 2-dim. dynamical motion
%Ref. PV-model in Brown & Hwang
% Brown, Robert Grover & Hwang, Patrick Y. C. (1997). Introduction to Random Signals
% and Applied Kalman Filtering. John Wiley & Sons, Inc., New York, 3rd edition.

%Q = Sv*[tau^3/3, tau^2/2; tau^2/2, tau];
%tau = .5;
Q = Sv*[tau, tau^2/2; tau^2/2, tau];
Q = blkdiag(Q,Q);

%Allocate memory to store all tracks over time. Each element of X is
%initially empty, but will store a matrix with 4 rows and T(p,k) columns,
%the j-th column of the matrix will store the 2-D position and velocity
%state variables of the j-th partition, of the p-th particle, at time k.
X = cell(num_targs,kmax);
truth = cell(num_targs,kmax);
Xinit = cell(num_targs,1);
Pinit = cell(num_targs,1);

%Allocate memory to store state covariance over time
P = cell(num_targs,kmax);

for z = 1:mc_runs
    
    z
    
    %Initialize each track state and error covariance
    for p = 1:num_targs
        Xinit{p} = zeros(4,1);
        
        Xinit{p}(1,:) = Px(p,1);
        Xinit{p}(2,:) = Vx(p,1);
        Xinit{p}(3,:) = Py(p,1);
        Xinit{p}(4,:) = Vy(p,1);
        
        Pinit{p} = diag([1 1e-2 1 1e-2]);
        
        
    end
    
    % Initialize each track classification state
    tt = [.2 .2 .2 .2 .2 .2 .2 .2 .2 .2; ...
        .2 .2 .2 .2 .2 .2 .2 .2 .2 .2; ...
        .2 .2 .2 .2 .2 .2 .2 .2 .2 .2; ...
        .2 .2 .2 .2 .2 .2 .2 .2 .2 .2; ...
        .2 .2 .2 .2 .2 .2 .2 .2 .2 .2];
    % tt = [.2 .2 .2 0 .2 .2 .2 .2 .2 .2; ...
    %     .2 .2 .2 1 .2 .2 .2 .2 .2 .2; ...
    %     .2 .2 .2 0 .2 .2 .2 .2 .2 .2; ...
    %     .2 .2 .2 0 .2 .2 .2 .2 .2 .2; ...
    %     .2 .2 .2 0 .2 .2 .2 .2 .2 .2];
    
    % setup a full screen figure
    if mc_runs == 1
        hfig = figure(1);
        set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    
    for k = 1:kmax
        
        % plot the location of the sensor
        if mc_runs == 1
            plot(s(1), s(3), '^r', 'LineWidth', 1, 'MarkerSize',20);
            hold on;
        end
        
        % update each kinematic track estimate
        for p = 1:num_targs
            
            % simulate actual state
            truth{p,k} = [Px(p,k), Vx(p,k), Py(p,k), Vy(p,k)]';
            
            % predict state estimate
            % special case for first time step
            if k==1
                X{p,k} = F*Xinit{p};% + blkdiag(A,A)*C*randn(4,1);
                %P{p,1} = diag([100 100]);
            else
                X{p,k} = F*X{p,k-1};% + blkdiag(A,A)*C*randn(4,1);
            end
            
%             % plot the mean of each Gaussian track
%             if mc_runs == 1
%                 plot(X{p,k}(1), X{p,k}(3), '*b', 'LineWidth', 1);
%                 
%                 text(X{p,k}(1) + 5, X{p,k}(3) + 5, num2str(p));
%             end
            
            % predict covariance
            if k ==1
                P{p,k} = Pinit{p};
            else
                
                % use larger process noise if the track is an ET
                %             if ets(p) == 1
                P{p,k} = F * P{p,k-1} * F' + Q;
                %             else
                %                 P{p,k} = F * P{p,k-1} * F';
                %             end
            end
            
            if periodic_update && mod(k,time_step) == 0
                %if mod(k,time_step) == 0
                
                [H, R] = genHandR(X{p,k},s);
                
                % update state estimate (assume all tracks are measured)
                dx = truth{p,k}(1) - s(1);
                dy = truth{p,k}(3) - s(3);
                range = sqrt(dx^2 + dy^2);
                %truthm = [truth{p,k}(1) truth{p,k}(3)]';
                %truthm = atan2(dy,dx);
                truthm = [range atan2(dy,dx)]';
                
                %noise = [randn * .1 randn * .1]';
                %noise = randn * .1;
                noise = [randn * 1 randn * .05]';
                
                m = truthm + noise;
                m = [m(1)*cos(m(2)) m(1)*sin(m(2))]';
                res = m - H * X{p,k};%[X{p,k}(1)-s(1) 0 X{p,k}(3)-s(3) 0]';
                S = H * P{p,k} * H' + R;
                K = P{p,k} * H' * S^-1;
                X{p,k} = X{p,k} + K * res;
                
                % update state estimate covariance (assume all tracks are measured)
                P{p,k} = (eye(4) - K * H) * P{p,k};
                
                %end
            end
            
        end
        
        if mc_runs == 1
            tt
        end
        
        % if only updating the target with highest ERR
        if ~periodic_update && mod(k,time_step) == 0
            
            % cost matrix for each of 3 classes
%             cost_mat = [0 30 10 1 1;
%                 1 0 10 1 1;
%                 1 30 0 1 1;
%                 1 30 10 0  1;
%                 1 30 10 1 0];

            cost_mat = [0 30 20 10 1;
                1 0 20 10 1;
                1 30 0 10 1;
                1 30 20 0  1;
                1 30 20 10 0];

            % SRM NOTE: the below is a sum of expected values.  Think of each
            % track having a multinomial distribution with p given by the class
            % PMF and N given by the number of looks at that track.  The expected
            % value over all the tracks is the sum of the individuals, and then the
            % EV over all time steps is the sum of the sums.
            
            % update the type I and type II cost matrices
            % start by estimating the number of targets in each class
            totals = zeros(size(tt,1),1);
            for p = 1:num_targs
                
                totals = totals + tt(:,p) * hgram(p, k + (z-1) * kmax);
                
            end
            
            %[cost_mat] = cost_mat_build([.1 .1 .5 .0 .3], kmax, kmax - k, totals);
            
            
            % go through each sensor and calculate the expected risk reduction
            % (ERR) for each track
            err = zeros(num_targs, 1);
            ttc = tt;
            for p = 1:num_targs
                
                % make covariance symmetric (remove numerical inaccuracy issues)使协方差对称
                P{p,k} = P{p,k} * P{p,k}';
                P{p,k} = sqrtm(P{p,k});
                
                [H, R] = genHandR(X{p,k},s);
                
                Y = P{p,k}^-1;
                Y = Y + H' * R^-1 * H;
                tempP = Y^-1;
                
                % calculate ERR
                err(p) = calcERR(ttc(:,p), cm, cost_mat, P{p,k}([1 3],[1 3]), tempP([1 3],[1 3]));
                
                % calculate info gain
                %err(p) = calcInfoGain(ttc(:,p), cm, cost_mat, P{p,k}([1 3],[1 3]), tempP([1 3],[1 3]));
                
                % calculate KL divergence
                %err(p) = calcKLdiv(ttc(:,p), cm, cost_mat, P{p,k}([1 3],[1 3]), tempP([1 3],[1 3]));
                
            end
            
            % save off err for analysis
            err_collection(:,k/time_step) = err;
            
            % determine the index of the max ERR (myopic) and draw an indication
            % that the sensor chose to look here
            if mc_runs == 1
                err
            end
            [err_sorted, II] = sort(err, 'descend');

%             % use this for sequential tasking (i.e. no optimization)
%             update_order = circshift(update_order,-1);
%             II = update_order;

            % update the top n tracks in terms of ERR
            n = 1;
            
            %for I = II(1:n)'
            for I = II(1:num_targs)'
                
                if mc_runs == 1
                    I
                end
                
                p = I(1);
                
                
                % check that the truth is within FOV of the target.  If not,
                % don't update.
                diff = X{p,k} - truth{p,k};
                if abs(diff(1)) > 500 || abs(diff(3)) > 500
                    disp(['WARNING: Target ' num2str(p) ' not in FOV, not updating track.']);
                    continue;
                end
                
                % update histogram
                hgram(p, k + (z-1) * kmax) = hgram(p, k + (z-1) * kmax) + 1;
                
                if mc_runs == 1
                    % draw sensor field of view
                    hcap1=patch(X{p,k}(1) + sx, X{p,k}(3) + sy, 'g');
                    %set(hcap1,'FaceAlpha',0.5 );
                    
                    % draw line from sensor to aimpoint
                    plot([s(1) X{p,k}(1)], [s(2) X{p,k}(3)], '-k');
                end
                
                % distribution sample to generate measurements from the given tracks
                measure = classDistributionSample(cm, class_truthc(p));
                
                % bayesian update the classification track
                [tt, pm] = updateClassTracks( cm, measure, tt, p);
                
                % go ahead and take a measurement on this target track and update track
                % via Kalman filter equations
                [H, R] = genHandR(X{p,k},s);
                
                % update state estimate (assume all tracks are measured)
                dx = truth{p,k}(1) - s(1);
                dy = truth{p,k}(3) - s(3);
                range = sqrt(dx^2 + dy^2);
                %truthm = [truth{p,k}(1) truth{p,k}(3)]';
                %truthm = atan2(dy,dx);
                truthm = [range atan2(dy,dx)]';
                
                %noise = [randn * .1 randn * .1]';
                %noise = randn * .1;
                noise = [randn * 1 randn * 3.0462e-04]';
                
                m = truthm + noise;
                m = [m(1)*cos(m(2)) m(1)*sin(m(2))]';
                res = m - H * X{p,k};%[X{p,k}(1)-s(1) 0 X{p,k}(3)-s(3) 0]';
                S = H * P{p,k} * H' + R;
                K = P{p,k} * H' * S^-1;
                X{p,k} = X{p,k} + K * res;
                
                % update state estimate covariance (assume all tracks are measured)
                P{p,k} = (eye(4) - K * H) * P{p,k};
                
                break;
                
            end
        end
        
        % update plots
        for p = 1:num_targs
            
            if mc_runs == 1
                
                % plot all updated state estimates and covariances 绘制所有更新的状态估计和协方差
                plot(X{p,k}(1), X{p,k}(3), 'd', 'MarkerSize', 10, 'LineWidth', 3, 'Color', [max(tt(2:3,p)) max(tt([1 4 5],p)) 0]);
                
                % plot the true target name associated with this track 绘制与此轨道关联的真实目标名称
                text(X{p,k}(1) + 60, X{p,k}(3) + 60, num2str(p));
                
                [xd yd] = covariance_ellipse(P{p,k}([1 3],[1 3]), X{p,k}(1), X{p,k}(3));
                plot(xd, yd, '-r');
                
                
                % plot truth
                for kk = 1:k
                    %kk = k;
                    
                    % plot based on class
                    if class_truthc(p) == 2 || class_truthc(p) == 3
                        
                        plot(truth{p,kk}(1), truth{p,kk}(3), '.r');
                        
                    else
                        
                        plot(truth{p,kk}(1), truth{p,kk}(3), '.g');
                        
                    end
                end
                
                %[xd yd] = covariance_ellipse(P{p,k}([1 3],[1 3]), X{p,k}(1), X{p,k}(3));
                %plot(xd, yd, '-r');
            end
            
            % calc se
            diff = (X{p,k} - truth{p,k});
            mse(p,k + (z-1) * kmax) = diff' * diff;
        end
        if k==kmax
        saveas(gcf,'1.jpg');
        end
        if k==1||k==20||k==40||k==60

            saveas(gcf,'test.jpg');
            1;
        end
        
        if mc_runs == 1
            axis([-4000 4000 -3000 5000]);
            %axis equal;
            
            % Pause to let it render tracks, reports, and schedules
            
%             pause(0.1);
            
            % Save frame to movie
            frames = frames + 1;
            Fr(frames) = getframe(gcf);
            if make_movie
                aviobj = addframe(aviobj, Fr(frames));
            end
            
            clf
        end
    end
    
end

% lastly plot a bar chart of the mse values
figure(2);
bar(sqrt(mean(mse')));
title('RMSE');

% and show histogram
figure(3);
bar(sum(hgram')./mc_runs);
% title('Histogram of Average # of measurements per target');
title('每个目标的测量次数平均值直方图');
% ylabel(['Average # of target measurements over ' num2str(kmax/time_step) ' schedule times and ' num2str(mc_runs) ' mc runs']);
ylabel(['' num2str(kmax/time_step) '个时间片和' num2str(mc_runs) '次蒙特卡洛循环下的目标测量次数平均值']);
xlabel('目标编号');
% xlabel('target id');

% and show ERR/info gain changes
figure(4);
plot(err_collection');
legend('1','2','3','4','5','6','7','8','9','10');
% xlabel('Planning time step (sim-time / 2)')
% title('Expected Risk Reduction Vs. Planning time step');
% ylabel('Expected Risk Reduction')
xlabel('规划时间步长(sim-time / 2)')
title('预期风险降低 Vs.计划时间步骤');
ylabel('预期风险降低')


if make_movie
    aviobj = close(aviobj);
end
