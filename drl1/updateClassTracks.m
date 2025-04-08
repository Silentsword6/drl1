% TMout = output track matrix (columns are tracks)
% pm = ouput probability of measurement
% CM = confusion matrix (rows are measurements)
% midx = row of CM
% TM = track matrix (columns are tracks)
% tidx = column of track matrix
function [ TMout, pm ] = updateClassTracks( CM, midx, TM, tidx)

% setup output size
TMout = TM;

% track being measured
t1 = TM(:,tidx);

% measurement
m = CM(midx,:);

% probability of measurement
pm = m * t1;

% new probability distribution for t1 (sequential Bayes)
t1n = m' .* t1 / pm;

% update viewed track and the sum of all track probabilities
TMout(:,tidx) = t1n;

end
