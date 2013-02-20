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

%% define 'args' structure, with default arguments
args = struct(...
    'nSubjects', 1, ...
    'priorType', 'loglinear', ...
    'nBoots', 5, ...
    'plotStraps', false, ...
    'plotSmears', true, ...
    'fitLikeFun', true, ...
    'fitRates', true, ...
    'fitPriorFun', false, ...
    'fitLikeSpeed', false, ...
    'plotPMFs', true, ...
    'pmfFitType', 'normal', ...
    'normUnits', 'linear', ...
    'nSEs', 2);

%% confirm that the input is parameter-value pairs
assert(mod(nargin, 2) == 0, ...
    'runSpeedDiscriminationModel:parseSpeedDiscriminationInput:paramValuePairs', ...
    'Input arguments should be in parameter-value pairs');

%% get the parameter values
% number of subjects
if any(strcmpi('nsubjects', varargin))
    args.nSubjects = varargin{find(strcmpi('nsubjects', varargin)) + 1};
end

% prior type
if any(strcmpi('priortype', varargin))
    args.priorType = varargin{find(strcmpi('priortype', varargin)) + 1};
end

% number of bootstraps
if any(strcmpi('nboots', varargin))
    args.nBoots = varargin{find(strcmpi('nboots', varargin)) + 1};
end

% flag to plot the bootstraps themselves
if any(strcmpi('plotstraps', varargin))
    args.plotStraps = varargin{find(strcmpi('plotstraps', varargin)) + 1};
end

% plot standard errors as smears
if any(strcmpi('plotsmears', varargin))
    args.plotSmears = varargin{find(strcmpi('plotsmears', varargin)) + 1};
end

% fit functional form to the contrast-dependence of the likelihood
if any(strcmpi('fitlikefun', varargin))
    args.fitLikeFun = varargin{find(strcmpi('fitlikefun', varargin)) + 1};
end

% include the maximum firing rate of the likelihood as a parameter to be fit
if any(strcmpi('fitrates', varargin))
    args.fitRates = varargin{find(strcmpi('fitrates', varargin)) + 1};
end

%%% NOTE: Although I have included logical statements throughout the code
%%% for fitting both a functional form to the prior and the
%%% speed-dependence of the likelihood function, my code does not actually
%%% do so.

% fit a functional form to the prior itself
if any(strcmpi('fitpriorfun', varargin))
    args.fitPriorFun = varargin{find(strcmpi('fitpriorfun', varargin)) + 1};
end

% fit speed-dependence of likelihood function
if any(strcmpi('fitlikespeed', varargin))
    args.fitLikeSpeed = varargin{find(strcmpi('fitlikespeed', varargin)) + 1};
end

% plot psychometric functions
if any(strcmpi('plotpmfs', varargin))
    args.plotPMFs = varargin{find(strcmpi('plotpmfs', varargin)) + 1};
end

% functional form of the fitted psychometric functions
if any(strcmpi('pmffittype', varargin))
    args.pmfFitType = varargin{find(strcmpi('pmffittype', varargin)) + 1};
end

% units used to normalize the prior distribution
if any(strcmpi('normunits', varargin))
    args.normUnits = varargin{find(strcmpi('normunits', varargin)) + 1};
end

% number of standard errors included in the plots
if any(strcmpi('nses', varargin))
    args.nSEs = varargin{find(strcmpi('nses', varargin)) + 1};
end