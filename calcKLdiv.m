% Date Created: 2/22/13
% Description:
% This is a function to calculate the expected risk reduction when
% looking at a given track, generalized to more than 2 classes.
% Also includes kinematic risk.
%
% Assumptions: the "interest" classification is always the first row
% of tt and cm, "non-interest" classification is always the second row
function [ kldiv ] = calcInfoGain( tt, cm, cost_mat, Pbefore, Pafter)

% FOV radius
rad = 250;

% NOTE: this scales the gain by P(Target of interest)
% TOI is either class 2 or 3, so sum those value
p_toi = tt(2) + tt(3);

% return the KL divergence
m = zeros(size(Pbefore,1),1); % same mean, only care about covariance reduction
kldiv = p_toi * kl_gaussian(m, m, Pbefore, Pafter);


end

% cm = [p(1|1) p(1|2); p(2|1) p(2|2)];