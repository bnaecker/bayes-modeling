function L = loglike_loglinear(prs, ds)
%
% FUNCTION L = loglike_loglinear(prs, ds)
%
% INPUT:    prs - guesses or updated guesses of the parameters of the
%               locally log-linear prior model
%
% OUTPUT:   L - negative log-likelihood of the parameters given the data
%
% loglike_loglinear returns the log-likelihood of the parameters in prs given
% the psychophysical choice data in ds.data.T, under the assumption that
% the prior is locally log-linear in the neighborhood of the likelihood
% function.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

% Indices of the current subject and bootstrap
si = ds.info.currentSubject;
bi = ds.info.currentBoot;

% Assign the current parameter guesses
% SIGMA:    the variance of the likelihood distribution
% a:        the local slope of the log-linear prior
sigVals = prs(1 : ds.info.nUniqueContrasts);
a = prs(ds.info.nUniqueContrasts + 1 : end);

% Arrange the slopes of the prior at each velocity using interp1
% this is literally rearranging, but without using indices
refSlope = interp1(ds.data(si).refVels, a, ...
    ds.data(si).refV, 'linear');
testSlope = interp1(ds.data(si).refVels, a, ...
    ds.data(si).testV, 'linear', 'extrap');

% Arrange the likelihood widths
refSig = sigVals(ds.data(si).refCInds);
testSig = sigVals(ds.data(si).testCInds);

% Estimate distribution for locally log-linear prior
Z = (ds.data(si).testV - ds.data(si).refV + refSlope .* refSig .^ 2 ...
    - testSlope .* testSig .^ 2) ./ sqrt(refSig .^ 2 + testSig .^ 2);

% Take the cumulative normal of this distribution, bounded away from zero
EPS = 1e-12;
rho = min(max(normcdf(Z), EPS), 1 - EPS);

% Compute log-likelihood
L = -sum(ds.data(si).T(:, bi) .* log(rho) + ...
    (1 - ds.data(si).T(:, bi)) .* log(1 - rho));
