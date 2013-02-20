function args = parseSpeedDiscrimInput(varargin)
%
% args = parseSpeedDiscrimInput(arglist)
%
% Helper function that parses the input, which should be in param-value
% pairs, to set all the flags. The default values are set here as well.
%
% INPUT:    varargin
%
% OUTPUT:   args - a structure containing the values for each of the
%                  possible inputs as fields
%
% (c) bnaecker@stanford.edu 17 Sep 2012

%% Define 'args' structure, with default arguments
args = struct(...
    'nSubjects', 1, ...
    'priorType', 'loglinear', ...
    'nBoots', 5, ...
    'plotStraps', false, ...
    'plotSmears', true, ...
    'fitLikeFun', true, ...
    'fitRates', true, ...
    'fitPriorFun', false, ...
    'fitG', false, ...
    'plotPMFs', true, ...
    'pmfFitType', 'normal', ...
    'normUnits', 'linear', ...
    'nSEs', 2);

%% Check that we've got param-value pairs
if mod(nargin, 2) ~= 0
    error('RunSpeedDiscriminationModel:parseSpeedDiscrimInput:paramValuePairs', ...
        'Inputs should be in parameter-value pairs');
else
    params = varargin(1:2:end);
    vals = varargin(2:2:end);
end

%% Check for each pair
% number of subjects
if any(strcmpi('nsubjects', params))
    args.nSubjects = vals{strcmpi('nsubjects', params)};
end

% prior type
if any(strcmpi('priortype', params))
    args.priorType = vals{strcmpi('priortype', params)};
end

% number of bootstraps
if any(strcmpi('nboots', params))
    args.nBoots = vals{strcmpi('nboots', params)};
end

% plot the bootstraps
if any(strcmpi('plotstraps', params))
    args.plotStraps = vals{strcmpi('plotstraps', params)};
end

% plot errors as smears
if any(strcmpi('plotsmears', params))
    args.plotSmears = vals{strcmpi('plotstraps', params)};
end

% fit likelihood function
if any(strcmpi('fitlikefun', params))
    args.fitLikeFun = vals{strcmpi('fitlikefun', params)};
end

% fit firing rates (max and min) of likelihood function
if any(strcmpi('fitrates', params))
    args.fitRates = vals{strcmpi('fitrates', params)};
end

% fit prior function (LATER FUNCTIONALITY)
if any(strcmpi('fitpriorfun', params))
    args.fitPriorFun = vals{strcmpi('fitpriorfun', params)};
end

% fit speed dependence of likelihood
% NOTE: Although I have included logical statements throughout the code to
% fit the speed-dependence of the likelihood function, because of the
% evidence suggesting there is no speed-dependence, I have not actually
% included the functionality. That is, though there is space in the code
% for one to do so, the current code does not.
if any(strcmpi('fitg', params))
    args.fitG = vals{strcmpi('fitg', params)};
end

% plot psychometric functions
if any(strcmpi('plotpmfs', params))
    args.plotPMFs = vals{strcmpi('plotpmfs', params)};
end

% functional form of psychometric function
if any(strcmpi('pmffittype', params))
    args.pmfFitType = vals{strcmpi('pmffittype', params)};
end

% units used to normalize the prior distribution
if any(strcmpi('normunits', params))
    args.normUnits = vals{strcmpi('normunits', params)};
end

% number of standard errors in plots
if any(strcmpi('nses', params))
    args.nSEs = vals{strcmpi('nses', params)};
end