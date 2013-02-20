function ds = RunSpeedDiscriminationModel(varargin)
%
% FUNCTION ds = RunSpeedDiscriminationModel(nSubjects, priorType)
%
% INPUT:    nSubjects = the number of subjects on which to run the model
%           priorType = string describing the type of the model to fit
%               'loglinear' or 'gaussian', describing the prior
%
% OUTPUT:   ds = structure with fitted parameters
%               Plots can be requested in setupData.m
%
% This is the top-level function to fit a Bayesian model to speed
% discrimination data. Please see the included README for more information.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

%% Parse input
args = parseSpeedDiscrimInput(varargin{:});
% if nargin == 0
%     priorType = 'loglinear';
%     nSubjects = 2;
% elseif nargin == 1 && isnumeric(varargin{1})
%     nSubjects = varargin{1};
%     priorType = 'loglinear';
%     warning('RunSpeedDiscriminationModel:undefinedPriorType', ...
%         'Setting prior type to "loglinear".');
% elseif nargin == 1 && ischar(varargin{1})
%     priorType = varargin{1};
%     nSubjects = 2;
%     warning('RunSpeedDiscriminationModel:undefinedNumberOfSubjects', ...
%         'Setting nSubjects to 2.');
% elseif nargin == 2
%     nSubjects = varargin{1};
%     priorType = varargin{2};
% elseif nargin > 2
%     nSubjects = varargin{1};
%     priorType = varargin{2};
%     warning('RunSpeedDiscriminationModel:tooManyInputArgs', ...
%         ['"RunSpeedDiscriminationModel" accepts up to two arguments. '...
%         '"nSubjects" defines the number of subjects on which to fit the model, '...
%         'and "priorType" defines the type of prior to fit, either ' ...
%         '"loglinear" or "gaussian". Ignoring extra inputs.']);
% end

%% Setup data structure
% Load data and organize into structure 'ds'
ds = setupData(args);

%% Collect data
% For each subject, pull out the relevant data and organize
ds = collectData(ds);

%% Setup optimization
% Setup the fitting parameters of the likelihood and prior distributions,
% and, if requested, fitting a functional form to the contrast-dependence
% of the likelihood function
ds = setupOptimization(ds);

%% Fit the model
% Using MATLAB's fmincon function, fit the model by maximimum likelihood.
% This is run on both the original data and bootstrapped datasets.
ds = fitModel(ds);

%% Calculate statistics of the fits
ds = getStatistics(ds);

%% Fit functional form to contrast-dependence of likelihood
ds = fitLikelihoodFun(ds);

%% Plot the fits
ds = plotDistFits(ds);

%% Plot the fitted psychometric functions
ds = plotPsychometricFunctions(ds);

%% Calculate and plot the subjects' thresholds
% ds = plotThresholds(ds);