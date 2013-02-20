function ds = reconstructPrior(ds)
%
% Helper function to reconstruct the prior, either from the fitted prior
% slopes, in the case of 'loglinear' model, or from the fitted gaussian
% parameters, in the case of a 'gaussian' prior.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

%% Get current velocities and define speed axis
if strcmp(ds.flags.normUnits, 'linear')
    vels = ds.data(ds.info.currentSubject).refVels;
    maxVel = max(vels) + .5 * max(vels);
    
    % Set linear normalization domain
    dx = .01;
    ds.params(ds.info.currentSubject).interpAx = min(vels):dx:(maxVel + dx);
elseif strcmp(ds.flags.normUnits, 'log')
    
    vels = ds.data(ds.info.currentSubject).refVels;
    maxVel = max(vels) + .1 * max(vels);
    
    % Set logarithmic normalization domain
    dx = .01;
    ds.params(ds.info.currentSubject).interpAx = ...
        ds.vTrans.transFun(dx:dx:(maxVel + dx));
else
    error('RunSpeedDiscriminationModel:fitModel:reconstructPrior', ...
        'Unsupported normalization units. Supported types are "log" and "linear"');
end

%% Normalize, depending on the model
if strcmp(ds.flags.priorType, 'loglinear')
    
    % Get the slopes returned by the optimization
    slopes = ds.params(ds.info.currentSubject).slopeHat(ds.info.currentBoot, :);
    
    % Interpolate slopes for normalization
    ySlopes = interp1(vels, slopes, ds.params(ds.info.currentSubject).interpAx,...
        'linear', 'extrap');
    
    % Sum and normalize the prior itself
    yCumulative = cumsum(ySlopes) * dx;
    prior = exp(-yCumulative);
    ds.params(ds.info.currentSubject).interpPrior(ds.info.currentBoot, :) = ...
        prior ./ (sum(prior) * dx);
    ds.params(ds.info.currentSubject).prior(ds.info.currentBoot, :) = ...
        interp1(ds.params(ds.info.currentSubject).interpAx, ...
        ds.params(ds.info.currentSubject).interpPrior(...
        ds.info.currentBoot, :), vels, 'linear');
    
else % Gaussian prior, assume zero mean
    pp = ...
        normpdf(ds.params(ds.info.currentSubject).interpAx, 0, ...
        ds.params(ds.info.currentSubject).gamma(ds.info.currentBoot));
    ds.params(ds.info.currentSubject).interpPrior(ds.info.currentBoot, :) = ...
        pp ./ (sum(pp) * dx);
    ds.params(ds.info.currentSubject).prior(ds.info.currentBoot, :) = ...
        interp1(ds.params(ds.info.currentSubject).interpAx, ...
        ds.params(ds.info.currentSubject).interpPrior(ds.info.currentBoot, :), vels);
end