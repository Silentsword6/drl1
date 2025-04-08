% Date Created: 1/17/13
% Description:
% This is a function to calculate the expected risk reduction when
% looking at a given track, generalized to more than 2 classes.
% Also includes kinematic risk.
%
function [ err ] = calcERR( tt, cm, cost_mat, Pbefore, Pafter)

% FOV radius
rad = 250;

% calculate the cumulative density over the FOV
% Guassian centered at zero with the covariance prior to a covar update
%pmiss_before = 1.0 - mvncdf([-1 -1 -1], [1 1 1], [0 0 0], Pbefore);
pmiss_before = 1.0 - mvncdf([-rad -rad], [rad rad], [0 0], Pbefore);

% and do the same with the covariance after a covar update
%pmiss_after = 1.0 - mvncdf([-1 -1 -1], [1 1 1], [0 0 0], Pafter);
pmiss_after = 1.0 - mvncdf([-rad -rad], [rad rad], [0 0], Pafter);

% calculate the current risk R
dim = size(cost_mat, 1);
Rvals = zeros(dim,1);
for idx = 1:dim
    Rvals(idx) = ...
        (cost_mat(:,idx)' * ones(dim,1) * tt(idx) * pmiss_before) * (dim-1) / dim + ... % pt. A
        cost_mat(idx,:) * tt * pmiss_before * 1 / dim + ...  % pt. B
        cost_mat(idx,:) * tt * (1 - pmiss_before); % pt. C

end

R = min(Rvals);

% calculate the expected risk R' after a measurement on this track
Rp = 0;
for m = 1:dim
    Rvals = zeros(dim,1);
    for idx = 1:dim
        Rvals(idx) =  ...
            (cost_mat(:,idx)' * ones(dim,1) * tt(idx) * cm(m,idx) * pmiss_after) * (dim-1) / dim + ... % pt. A
            cost_mat(idx,:) * (tt .* cm(m,:)')  * pmiss_after * 1 / dim + ... % pt. B
            cost_mat(idx,:) * (tt .* cm(m,:)')  * (1 - pmiss_after); % pt. C
            
    end
    Rp = Rp + min(Rvals);
end

% return the expected risk reduction
err = R - Rp;

end

% cm = [p(1|1) p(1|2); p(2|1) p(2|2)];
