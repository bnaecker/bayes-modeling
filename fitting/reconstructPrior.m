function ds = reconstructPrior(ds)
%
% FUNCTION ds = reconstructPrior(ds)
%
% Helper function to reconstruct the prior distribution itself, either from
% the fitted prior slopes in the case of the loglinear model, or from the
% fitted variances in the case of a gaussian prior.
% 
% (c) bnaecker@stanford.edu 11 Nov 2012

%% counters for the current subject and the current bootstrap iteration
si = ds.info.currentSubject;
bi = ds.info.currentBoot;

%% get current velocities and define the correct speed axis
if strcmp(ds.flags.normUnits, 'linear')
    vels = ds.data(si).refVels;
    maxVel = max(vels) + 0.5 * max(vels);
    
    % set linear normalization units
    dx = 0.01;
    ds.params(si).interpAx = dx:dx:(maxVel + dx);
    
else
    vels = ds.data(si).refVels;
    maxVel = max(vels) + 0.1 * max(vels);
    
    % set logarithmic normalization units
    dx = 0.01;
    ds.params(si).interpAx = ...
        ds.velTrans.transFun(dx : dx : (maxVel + dx));
end

%% normalize, by model
switch ds.flags.priorType
case 'loglinear'
    % get the prior slopes returned by the minimization procedure
    slopes = ds.params(si).slopeHat(:, bi);
    
    % interpolate the slopes for the normalization
    ySlopes = [interp1(vels, slopes, ds.params(si).interpAx(ds.params(si).interpAx < vels(1)), ...
		'linear', slopes(1)) ...
	interp1(vels, slopes, ds.params(si).interpAx(ds.params(si).interpAx >= vels(1) & ...
		ds.params(si).interpAx <= vels(end)), 'linear') interp1(vels, slopes, ...
		ds.params(si).interpAx(ds.params(si).interpAx > vels(end)), 'linear', slopes(end))];
    
    % integrate the slopes and normalize the resulting distribution itself
    yCumulative = cumsum(ySlopes) * dx;
    prior = exp(-yCumulative);
    ds.params(si).interpPrior(:, bi) = ...
        prior ./ (sum(prior) * dx);
    
    % interp1 is used to pick the slopes from the normalized distribution
    % at the values of the reference velocities
    ds.params(si).prior(:, bi) = ...
        interp1(ds.params(si).interpAx, ...
        ds.params(si).interpPrior(:, bi), ...
        vels, 'linear');
    
case 'gaussian'
    % compute the normal distribution function with fitted variance
    pp = normpdf(ds.params(si).interpAx, 0, ...
        ds.params(si).gamma(bi) ^ 2);
    
    % normalize
    ds.params(si).interpPrior(:, bi) = ...
        pp ./ (sum(pp) * dx);
    
    % use interp1 to pick out the values of the distribution at the
    % reference velocities
    ds.params(si).prior(:, bi) = ...
        interp1(ds.params(si).interpAx, ...
        ds.params(si).interpPrior(:, bi), vels);

case 'mixture'
	% compute mixture distribution, individual components summed 
	pp = sum(normpdf(ds.params(si).interpAx' * ...
		ones(1,ds.mixInfo.nComponents), ...
		ones(length(ds.params(si).interpAx, 1), 1) * ...
		ds.params(si).mixMeans, ...
		ones(length(ds.params(si).interpAx, 1), 1) * ...
		ds.params(si).mixMeans), 2);

	% normalize
	ds.params(si).interpPrior(:, bi) = ...
		pp ./ (sum(pp) * dx);

	% use interp1 to pick out the values of the distribution at the reference velocities
	ds.params(si).prior(:, bi) = ...
		interp1(ds.params(si).interpAx, ...
		ds.params(si).interpPrior(:, bi), vels);

otherwise 
    % optional functionality
end
