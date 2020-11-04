function args = parseSpeedDiscriminationInput(varargin)
%
% FUNCTION args = parseSpeedDiscriminationInput(varargin)
%
% Helper function that parses the requested parameter-value pairs, and sets
% up the data structure.
%
% INPUT:    varargin - parameter-value pairs
%
% OUTPUT:   args - structure containing fields with requested or default
%               parameter values.
%
% (c) bnaecker@stanford.edu 11 Nov 2012
% 2017-01-11 - update to use inputParser

parser = inputParser();
arguments = struct( ...
	'nSubjects', ...	% Number of subjects used
		{{1, @(x) x > 0}}, ...
	'priorType', ...	% Functional form of prior to be fitted
		{{'loglinear', @(x) ischar(x) && (any(strcmp(x, {'loglinear', 'gaussian', 'mixture'})))}}, ...
	'nBoots', ...		% Number of bootstrap iterations on fitting procedure
		{{5, @(x) x > 0}}, ...
	'plotStraps', ...	% Plot each individual bootstrap, where appropriate
		{{false, @(x) islogical(x)}}, ...
	'plotSmears', ...	% Plot polygons showing errors, instead of individual bootstraps
		{{true, @(x) islogical(x)}}, ...
	'fitLikeFun', ...	% Fit functional form to contrast-dependence of likelihood function
		{{true, @(x) islogical(x)}}, ...
	'fitRates', ...		% Include maximum firing rate in functional form of likelihood
		{{true, @(x) islogical(x)}}, ...
	'fitLikeSpeed', ...	% Fit speed-dependence of likelihoods
		{{false, @(x) islogical(x)}}, ...
	'plotPMFs', ...		% Plot psychometric functions
		{{true, @(x) islogical(x)}}, ...
	'pmfFitType', ...	% Type of function used to fit PMFs
		{{'normal', @(x) ischar(x) && (any(strcmp(x, {'normal', 'weibull'})))}}, ...
	'normUnits', ...	% Units used when normalizing the prior distribution
		{{'linear', @(x) ischar(x) && (any(strcmp(x, {'linear', 'log'})))}}, ...
	'nSEs', ...			% Number of standard errors to plot
		{{2, @(x) x >= 0.}}, ...
	'nPmfSamples', ...	% Number of samples drawn at each velocity when fitting PMFs
		{{100, @(x) x > 0}}, ...
	'nComponents', ...	% Number of mixture components
		{{5, @(x) x >= 1}});

argnames = fieldnames(arguments);
for i = 1:length(argnames)
	parser.addOptional(argnames{i}, ...
		arguments.(argnames{i}){1}, ... % default value
		arguments.(argnames{i}){2}); 	% validator function
end

% Parse
parser.parse(varargin{:});
args = parser.Results;
