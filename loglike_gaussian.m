function L = loglike_gaussian(prs, ds)
%
% FUNCTION L = loglike_gaussian(prs, ds)
%
% INPUT:    prs = guesses or updated guesses of the parameters of the full
%               gaussian prior distribution
%
% OUTPUT:   L = log-likelihood of the parameters given the data
%
% logli_gaussian returns the log-likelihood of the parameters in prs given
% the psychophysical choice data in data.data.T, under the assumption that
% the prior is a gaussian distribution
%
% (c) bnaecker@stanford.edu 11 Feb 2012

% Indices of the current subject and bootstrap
si = ds.info.currentSubject;
bi = ds.info.currentBoot;

% Assign the current parameter guesses
% SIGMA:    the variance of the likelihood distribution
% GAMMA:    the variance of the prior distribution
sigVals = prs(1 : ds.info.nUniqueContrasts);
gamma = prs(end);

% Arrange the likelihood variances
refSig = sigVals(ds.data(si).refCInds);
testSig = sigVals(ds.data(si).testCInds);

% Compute alpha = (gamma^2 / (gamma^2 + sigma^2) at each test velocity
alpha1 = gamma .^ 2 ./ (gamma .^ 2 + refSig .^ 2);
alpha2 = gamma .^ 2 ./ (gamma .^ 2 + testSig .^ 2);

% Compute the estimate distribution for the gaussian prior
Z = (alpha2 .* ds.data(si).testV - alpha1 .* ds.data(si).refV) ./ ...
    sqrt( alpha1 .^ 2 .* (refSig .^ 2) + alpha2 .^ 2 .* (testSig .^ 2) );

% Compute estimate distribution (PMF)
EPS = 1e-12;
rho = min(max(normcdf(Z), EPS), 1 - EPS);

% Compute log-likelihood
L = sum(ds.data(si).T(:, bi) .* log(rho) + ...
    (1 - ds.data(si).T(:, bi)) .* log(1 - rho));