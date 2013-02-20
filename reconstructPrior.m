function ds = reconstructPrior(ds)
%
% FUNCTION ds = reconstructPrior(ds)
%
% Helper function to reconstruct the prior distribution itself, either from
% the fitted prior slopes in the case of the loglinear model, or from the
% fitted variances in the case of a gaussian prior.
% 
% (c) bnaecker@stanford.edu 11 Nov 2012

%% get current velocities and define the correct speed axis
if strcmp(ds.flags.normUnits, 'linear')
    vels = ds.data(ds.info.currentSubject).refVels;
    maxVel = max(vels) + 0.5 * max(vels);
    
    % set linear normalization units
    dx = 0.01;
    ds.params(ds.info.currentSubject).interpAx = min(vels):dx:(maxVel + dx);
    
else
    vels = ds.data(ds.info.currentSubject).refVels;
    maxVel = max(vels) + 0.1 * max(vels);
    
    % set logarithmic normalization units
    dx = 0.01;
    ds.params(ds.info.currentSubject).interpAx = ...
        ds.velTrans.transFun(dx : dx : (maxVel + dx));
end

%% normalize, by model
if strcmp(ds.flags.priorType, 'loglinear');
    % get the prior slopes returned by the minimization procedure
    slopes = ds.params(ds.info.currentSubject).slopeHat(:, ds.info.currentBoot);
    
    % interpolate the slopes for the normalization
    ySlopes = interp1(vels, slopes, ds.params(ds.info.currentSubject).interpAx, ...
        'linear', 'extrap');
    
    % integrate the slopes and normalize the resulting distribution itself
    yCumulative = cumsum(ySlopes) * dx;
    prior = exp(-yCumulative);
    ds.params(ds.info.currentSubject).interpPrior(:, ds.info.currentBoot) = ...
        prior ./ (sum(prior) * dx);
    
    % interp1 is used to pick the slopes from the normalized distribution
    % at the values of the reference velocities
    ds.params(ds.info.currentSubject).prior(:, ds.info.currentBoot) = ...
        interp1(ds.params(ds.info.currentSubject).interpAx, ...
        ds.params(ds.info.currentSubject).interpPrior(:, ds.info.currentBoot), ...
        vels, 'linear');
    
elseif strcmp(ds.flags.priorType, 'gaussian')
    % compute the normal distribution function with fitted variance
    pp = normpdf(ds.params(ds.info.currentSubject).interpAx, 0, ...
        ds.params(ds.info.currentSubject).gamma(ds.info.currentBoot) ^ 2);
    
    % normalize
    ds.params(ds.info.currentSubject).interpPrior(:, ds.info.currentBoot) = ...
        pp ./ (sum(pp) * dx);
    
    % use interp1 to pick out the values of the distribution at the
    % reference velocities
    ds.params(ds.info.currentSubject).prior(:, ds.info.currentBoot) = ...
        interp1(ds.params(ds.info.currentSubject).interpAx, ...
        ds.params(ds.info.currentSubject).interpPrior(:, ds.info.currentBoot), vels);
    
else 
    % optional functionality
end