function L = loglike_mixture(prs, ds)
%
% FUNCTION L = loglike_mixture(prs, ds)
%
% INPUT:    prs = guesses or updated guesses of the parameters of the full
%               gaussian prior distribution
%
% OUTPUT:   L = negative log-likelihood of the parameters given the data
%
% logli_mixture returns the log-likelihood of the parameters in prs given
% the psychophysical choice data in ds.data.T, under the assumption that
% the prior is a mixture of gaussians distribution
%
% (c) bnaecker@stanford.edu 03 Apr 2013

%% indices of the current subject and bootstrap
si = ds.info.currentSubject;
bi = ds.info.currentBoot;

%% assign the current parameter guesses
% variance of the likelihood distribution
sigVals = prs(1 : ds.info.nUniqueContrasts);

% mixture component means
mu = prs(ds.info.nUniqueContrasts + 1 : ...
	ds.info.nUniqueContrasts + ds.mixInfo.nComponents);

% variances of the mixture components
gamma = prs(ds.info.nUniqueContrasts + ds.mixInfo.nComponents + 1 : ...
	ds.info.nUniqueContrasts + 2 * ds.mixInfo.nComponents);

% mixing weights for each component
w = prs(ds.info.nUniqueContrasts + 2 * ds.mixInfo.nComponents + 1 : end);

%% arrange the likelihood variances
refSig = sigVals(ds.data(si).refCInds);
testSig = sigVals(ds.data(si).testCInds);

%% compute parameters of estimate distribution
% this includes the shrinkage factors, modified means (mu tilde in the paper), and
% modified mixing weights (w tilde). these must be computed for both the 
% reference and test stimuli, on each trial

% preallocate
refAlpha = zeros(ds.mixInfo.nComponents, ds.info.nTrials);
testAlpha = refAlpha;
refMu = refAlpha;
testMu = refAlpha;
refW = refAlpha;
testW = refAlpha;

% loop over components, computing the required values
for ci = 1:ds.mixInfo.nComponents
	refAlpha(ci, :) = (gamma(ci)^2 ./ (refSig^2 + gamma(ci)^2));
	testAlpha(ci, :) = (gamma(ci)^2 ./ (testSig^2 + gamma(ci)^2));
	refMu(ci, :) = (refSig^2 / (refSig^2 + gamma(ci)^2)) * mu(ci);
	testMu(ci, :) = (testSig^2 / (testSig^2 + gamma(ci)^2)) * mu(ci);
	refV = (w(ci) / sqrt(gamma(ci)^2 + refSig^2)) * normpdf( (m - refMu(ci, :)) / sqrt(gamma(ci)^2 + refSig^2));
	testV = (w(ci) / sqrt(gamma(ci)^2 + testSig^2)) * normpdf( (m - testMu(ci, :)) / sqrt(gamma(ci)^2 + testSig^2));

	%wRef(ci, :) = (w(ci) .* sqrt(gamma(ci)^2 + refSig.^2)) .* ...
		%normpdf(
	%% PROBLEM! WE DON'T HAVE M!
end

% compute modified means, for reference and test velocities, on each trial
for ci = 1:ds
end
