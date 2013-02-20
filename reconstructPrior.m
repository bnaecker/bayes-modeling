function var = reconstructPrior(var)
%
% helper function var = reconstructPrior(var) takes the returned prior
% slopes from fmincon and generates the actual prior probability
% distribution
%
% (c) benjamin naecker UT Austin 8 Jun 2011 benjamin.naecker@gmail.com

if strcmp(var.flags.normUnits, 'linear')
    
    % Get current velocities and slopes
    vels = var.vTrans.iTransFun(var.data(var.info.currentSubject).refVels);
    slopes = var.params(var.info.currentSubject).slopeHat(var.info.currentBoot, :);
    
    % Set normalization domain
    dx = .01;
    xx = min(vels):dx:(max(vels) + dx);
    
elseif strcmp(var.flags.normUnits, 'log')
    
    % Get current velocities and slopes
    vels = var.data(var.info.currentSubject).refVels;
    slopes = var.params(var.info.currentSubject).slopeHat(var.info.currentBoot, :);
    
    % Set normalization domain
    dx = .01;
    xx = var.vTrans.iTransFun(var.vTrans.v0):dx:(max(vels) + dx);
    
else
    error('Unsupported normalization for reconstructing the prior from slopes');
end

% Interpolate slopes for normalization
ySlopes = interp1(vels, slopes, xx, 'cubic', 'extrap');

% Sum and normalize prior
yCum = cumsum(ySlopes)*dx;
prior = exp(-yCum);
var.params(var.info.currentSubject).interpPrior(var.info.currentBoot, :) = ...
    prior ./ sum(prior)/dx;
var.params(var.info.currentSubject).prior(var.info.currentBoot, :) = ...
    interp1(xx, var.params(var.info.currentSubject).interpPrior(var.info.currentBoot, :), ...
    vels, 'linear');