function ds = setupData(args)
%
% FUNCTION ds = setupData(args)
%
% INPUT:    args = the structure output of parseSpeedDiscrimInput.m
%
% OUTPUT:   ds = data structure containing everythang
%
% Loads data files from disk containing psychophysical data collected by AA
% Stocker and EP Simoncelli. This function assumes that the data is
% contained in a folder within this directory, called 'Data', and that the
% subjets' *.mat files are titled 's1.mat', etc.
%
% Also described here are the log-transform of the velocity axis and all
% flags.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

%% Initialize data structure
ds = struct(...
    'data', {[]}, ...       % Contains all the subjects' data
    'flags', {[]}, ...      % Flags for fitting, plotting, etc
    'hOpt', {[]}, ...       % Optimization setup for fitting likelihood
    'handles', {[]}, ...    % Handles to all plots
    'info', {[]}, ...       % General information about the model
    'llOpt', {[]}, ...      % Optimization setup for the log-likelihood
    'params', {[]}, ...     % Parameters returned by the optimization
    'stats', {[]}, ...      % Statistics of resampled fits
    'pmfs', {[]}, ...       % Psychometric functions (model-predicted and/or Weibull fits)
    'vTrans', {[]});        % Info for the transformation of the velocity axis

%% Setup the file structure and load data
fprintf('\nFinding files and loading...');
ds.info.baseDir = pwd;
ds.info.dataDir = fullfile(ds.info.baseDir, 'Data');
if args.nSubjects == 0
    ds.info.nSubjects = 1;
    if strcmp(args.priorType, 'gaussian')
        d = load(fullfile(ds.info.dataDir, 's3.mat'));
        field = 's3';
    elseif strcmp(args.priorType, 'loglinear')
        d = load(fullfile(ds.info.dataDir, 's4.mat'));
        field = 's4';
    else
        error('RunSpeedDiscriminationModel:setupData:unknownPriorType', ...
            'Supported prior types are "loglinear" and "gaussian"');
    end
    ds.data(1).rawData = d.(field);
    ds.data.signse = [d.sig1 d.sig2];
    ds.info.isSimData = true;
else
    ds.info.nSubjects = args.nSubjects;
    ds.info.isSimData = false;
    for si = 1:ds.info.nSubjects
        ds.info.dataFiles{si} = ['s' num2str(si) '.mat'];
        d = load(fullfile(ds.info.dataDir, ds.info.dataFiles{si}));
        fieldName = ['s' num2str(si)];
        ds.data(si).rawData = d.(fieldName)';
    end
end
fprintf('done.\nSetting up parameters and flags...');

%% Organize column indices
% These are specific to the data collected from AA Stocker and EP
% Simoncelli. Tailor to your data sets.
ds.info.refVCol = 1;        % Speed of reference stimulus
ds.info.refCCol = 2;        % Contrast of reference stimulus
ds.info.testVCol = 3;       % Speed of test stimulus
ds.info.testCCol = 4;       % Contrast of test stimulus
ds.info.stairCol = 5;       % Indicates which of two interleaved staircases
ds.info.dirCol = 6;         % Direction of stimuli
ds.info.placeCol = 7;       % Location of reference stimulus (1 == right)
ds.info.rtCol = 8;          % Reaction time
ds.info.choiceCol = 9;      % Subject's choice (1 == test seen faster)
ds.info.adaptCol = 10;      % Step direction of staircase
ds.info.correctCol = 11;    % Correct choice

%% Setup log-transform of velocity axis
% Stocker and Simoncelli assume that representation of the likelihood
% function is normal on a logarithmic speed axis. This is borne out by
% electrophysiological evidence (See their Methods for discussion).
ds.vTrans.v0 = .3;        % v0
ds.vTrans.transFun = ...  % Tranformation to log-axis
    @(v) ( log(1 + (v ./ ds.vTrans.v0)) );
ds.vTrans.iTransFun = ... % Inverse log-transform
    @(v) ( expm1(v) .* ds.vTrans.v0 );

%% Flags
% The default values and explanations for these are set in
% parseSpeedDiscrimInput.m

% collect everything but the number of subjects
ds.flags = args;
ds.flags = rmfield(ds.flags, {'nSubjects', 'nBoots', 'nSEs'});
ds.info.nBoots = args.nBoots; ds.info.nSEs = args.nSEs;

%% Define g(v), the speed-dependence of the likelihood function
if ~ds.flags.fitG
    if ds.info.nSubjects == 0
        ds.data(1).G = 1;
    else
        ds.data(1).G = .25;
        if ds.info.nSubjects > 1
            ds.data(2).G = .2;     % These were just read off figure 4
        end
    end
end

%% Define percent correct for plotting thresholds
ds.info.pCor = .82;

%% Organize
ds = orderDataStruct(ds);
fprintf('done.');