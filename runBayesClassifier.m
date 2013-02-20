% Run inverse Bayesian modeling method on data from Stocker & Simoncelli
% Nat Neurosci 2006 paper.
%
% This version is more streamlined, with subfunctions and better data
% structures and comments
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

%% Setup
% Load data and organize into 'var' structure
var = setupVar;

%% Organize data
% For each subject, pull out the relevant data and organize it for the
% log-likelihood function
var = organizeData(var);

%% Setup otpimization
% Setup the parameters of the fitting of the likelihood and the prior, and,
% if requested in var.flags, the functional form of the likelihood h(c)
var = setupOptimization(var);

%% Run optimization
% Run fmincon on negative log-likelihood function for original and
% bootstrapped data sets, returning the likelihood h(c) and the prior that
% the subject used
var = fitData(var);
if isfield(var, 'error')
    return
end

%% Calculate statistics of bootstrapped samples
% Calculate +/- 2 SEs from mean of likelihood and prior
var = getStats(var);

%% Fit functional form of h(c)
if var.flags.fitLikeFun
    var = fitLikeFun(var);
end

%% Plot fits to likelihood and prior, plus h(c) and individual bootstraps,
% if requested in var.flags
var = plotFits(var);

%% Plot subjects' psychometric functions, if request in var.flags
% This needs to run regardless of whether or not we plot, since the
% proportions calculated here and used in the PMFs are also used in
% calculating the thresholds below
% var = plotPMFs(var);

%% Calculate and plot subjects' thresholds
% both absolute (\Delta V) and relative thresholds (\Delta V / V)
% var = plotThresholds(var);