function L = logli_loglinear(prs, ds)
%
% FUNCTION L = logli_loglinear(prs, ds)
%
% INPUT:    prs = guesses or updated guesses of the parameters of the
%               locally log-linear prior model
%
%           ds = the main data structure
%
% OUTPUT:   L = log-likelihood of the parameters given the data
%
% logli_loglinear returns the log-likelihood of the parameters in prs given
% the psychophysical choice data in data.data.T, under the assumption that
% the prior is locally log-linear in the neighborhood of the likelihood
% function.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

% Indices of the current subject and bootstrap
subject = ds.info.currentSubject;
boot = ds.info.currentBoot;

% Assign the current parameter guesses
% SIGMA:    the variance of the likelihood distribution
% a:        the local slope of the log-linear prior
sigVals = prs(1 : ds.info.nUniqueContrasts);
a = prs(ds.info.nUniqueContrasts + 1 : end);

% Arrange the slopes of the prior at each velocity using interp1
refSlope = interp1(ds.data(subject).refVels, a, ...
    ds.data(subject).refV, 'linear');
testSlope = interp1(ds.data(subject).refVels, a, ...
    ds.data(subject).testV, 'linear', 'extrap');

% Arrange the likelihood widths
refSig = sigVals(ds.data(subject).refCInds);
testSig = sigVals(ds.data(subject).testCInds);

% Estimate distribution for locally log-linear prior
Z = (ds.data(subject).testV - ds.data(subject).refV + refSlope .* refSig .^ 2 ...
    - testSlope .* testSig .^ 2) ./ sqrt(refSig .^ 2 + testSig .^ 2);

% Take the cumulative normal of this distribution, bounded away from zero
EPS = 1e-12;
rho = min(max(normcdf(Z), EPS), 1 - EPS);

% Compute log-likelihood
L = sum(ds.data(subject).T(:, boot) .* log(rho) + ...
    (1 - ds.data(subject).T(:, boot)) .* log(1 - rho));