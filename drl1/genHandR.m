% generate H and R matrices given the position of the track and sensor
function [H, R] = genHandR(x, s)

%sensor azimuth and range measurement variance
variance_az = 3.0462e-04;
variance_rng = 1;

% pull out resulting x, y for simplicity
lx = x(1) - s(1);
ly = x(3) - s(3);

% calculate the range
range = sqrt(lx^2 + ly^2);

model = 4;

% 1 : bearing only
% 2 : euclidean position measure
% 3 : bearing and range
% 4 : bearing and range (linearize about noise)
% 5 : bearing only (linearize about noise)
if model == 1
    
    H = zeros(1,4);
    
    % calculate Jacobian
    H(1,1) = -ly / (ly*ly + lx*lx);
    H(1,3) = lx / (ly*ly + lx*lx);
    
    R = zeros(1,1);
    R(1,1) = variance_az;
    
elseif model == 2
    
    H = zeros(2,4);
    H(1,1) = 1;
    H(2,3) = 1;
    
    R = zeros(2,2);
    R(1,1) = variance_rng;
    R(2,2) = variance_rng;
    
elseif model == 3
    
    H = zeros(2,4);
    
    % calculate Jacobian
    H(1,1) = -ly / (ly*ly + lx*lx);
    H(1,3) = lx / (ly*ly + lx*lx);
    H(2,1) = lx / range;
    H(2,3) = ly / range;
    
    R = zeros(2,2);
    R(1,1) = variance_az;
    R(2,2) = variance_rng;
    
elseif model == 4
    
    heading = atan2(ly, lx);
    
    if heading > 2 * pi
        heading = heading - 2 * pi;
    end
    
    if heading < 0
        heading = heading + 2 * pi;
    end
    
    s_head = sin(heading);
    c_head = cos(heading);
    
    % setup the measurement covariance (R) and observation (H) matrices for
    % the sensor
    % NOTE: currently using a standard conversion to linearize the measurement
    % for R
    H = zeros(2,4);
    H(1,1) = 1;
    H(2,3) = 1;
    
    R = zeros(2,2);
    R11 = range^2 * variance_az * s_head^2 + variance_rng * c_head^2;
    R22 = range^2 * variance_az * c_head^2 + variance_rng * s_head^2;
    R12 = (variance_rng - range^2 * variance_az) * s_head * c_head;
    
    R(1,1) = R11; R(1,2) = R12; R(2,1) = R12; R(2,2) = R22;
    
    
end
