function ds = runSpeedDiscriminationModel(varargin)
%
% FUNCTION ds = runSpeedDiscriminationModel(varargin)
%
% This is the top-level function used to fit a Bayesian ideal observer
% model to discrimination data. The function takes arguments in
% parameter-value pairs, all of which are described in detail in the
% included README file.
%
% INPUT:    arguments in parameter-value pairs. see README for more info.
%
% OUTPUT:   ds - generic data structure including all raw data, model fits,
%               handle objects, and fit objects.
%
% (c) bnaecker@stanford.edu 11 Nov 2012


%% parse input
args = parseSpeedDiscriminationInput(varargin{:});

%% setup data structure
ds = setupDataStruct(args);

%% collect raw data
ds = collectRawData(ds);

%% setup optimization parameters
ds = setupOptimization(ds);

%% fit Bayesian ideal observer model
ds = fitModel(ds);

%% compute some statistics about the fits
ds = getFitStatistics(ds);

%% fit functional form to the contrast-dependence of likelihood function
ds = fitLikelihoodFunction(ds);

%% plot the fitted prior and likelihood widths
ds = plotPriorFit(ds);

%% plot the actual distributions with parameters fitted by the model
% ds = plotFittedDistributions(ds);

%% plot fitted psychometric functions
% ds = plotPsychometricFunctions(ds);