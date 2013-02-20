function var = setupVar
%
% setupData loads the files from disk containing the data for both subjects
% from the Stocker & Simoncelli paper, and organizes it and some basic
% indexing into the main struct 'var'
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

%% Initialize 'var' struct
var = struct( ...
    'data', [], ...
    'flags', [], ...
    'hOpt', [], ...
    'handles', [], ...
    'info', [], ...
    'llOpt', [], ...
    'params', [], ...
    'stats', [], ...
    'vTrans', []);

%% Load data files from disk and organize
% Flag which datasets to use
var.flags.useStockerData = 1;   % Use data from Stocker/Simoncelli paper (subjects 1 & 2)
var.flags.usePillowData = 0;    % Use data from Pillow lab

% Data from Stocker & Simoncelli 2006
if var.flags.useStockerData
    var.info.stockerBaseDir = '/Volumes/LKCLAB/Users/BenN/Code/BayesianIdealObserver/StockerData/';
    var.info.stockerBaseDir = '/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/Data/Stocker';
    var.info.stockerDataFiles{1} = 's1.mat';
    var.info.stockerDataFiles{2} = 's2.mat';
    load(fullfile(var.info.stockerBaseDir, var.info.stockerDataFiles{1}));
    load(fullfile(var.info.stockerBaseDir, var.info.stockerDataFiles{2}));
    var.data(1).rawData = s1';
    var.data(2).rawData = s2';
    nStockerSubjects = 2;
else
    nStockerSubjects = 0;
end

% Data from Pillow lab psychophysics
if var.flags.usePillowData
    % Need to fineagle it into the same format as Stocker data
    var.info.pillowBaseDir = '/Volumes/LKCLAB/Psychophysics/SpeedDiscrim/';
    ext = 'Data';
    priorType = 'high';
    subjectList = {'BNN'};
    for subject = 1:length(subjectList)
        var.info.pillowDataFiles{subject} = subjectList{subject};
        d = dir(fullfile(var.info.pillowBaseDir, ext, priorType, var.info.pillowDataFiles{subject}));
        load(fullfile(var.info.pillowBaseDir, ext, priorType, var.info.pillowDataFiles{subject}, ...
            d(end).name));
        var.data(subject + nStockerSubjects).rawData = ...
            [data(:, 3) data(:, 4) data(:, 8) data(:, 5) nan(size(data, 1), 1) data(:, 7) ...
            data(:, 6) nan(size(data, 1), 1) data(:, 9) nan(size(data, 1), 1) data(:, 10)]; %#ok<NODEF>
    end
    nPillowSubjects = 1;
else
    nPillowSubjects = 0;
end

%% Organize data
var.info.nSubjects = nStockerSubjects + nPillowSubjects;
               
%% Organize column indices
var.info.refVCol = 1;        % speed of reference
var.info.refCCol = 2;        % contrast of reference
var.info.testVCol = 3;       % speed of test
var.info.testCCol = 4;       % contrast of test
var.info.stairCol = 5;       % label of staircase (unused)
var.info.dirCol = 6;         % direction of gratings (unused)
var.info.placeCol = 7;       % location of reference grating (1 = right aperture)
var.info.rtCol = 8;          % reaction time (unused)
var.info.choiceCol = 9;      % subject's choice (1 = test seen faster)
var.info.adaptCol = 10;      % step direction of stair (unused)
var.info.correctCol = 11;    % correct choice

%% Setup log-transform function of speeds
var.vTrans.v0 = .3;
var.vTrans.transFun = @(v) ( log(1 + (v ./ var.vTrans.v0)) );
var.vTrans.iTransFun = @(v) ( expm1(v) .* var.vTrans.v0 );

%% Number of bootstrap iterations
var.info.nBoots = 5;

%% Number of standard errors in plots
var.info.nSEs = 2;

%% Standard deviation threshold for throwing out resampled fits
var.info.SDThresh = 2.8;

%% Flags
var.flags.plotStraps = 1;       % Plot each bootstrap iteration
if var.info.nBoots == 1
    var.flags.plotSmears = 0;
else
    var.flags.plotSmears = 1;   % Plot SEs on output data
end
var.flags.fitLikeFun = 1;       % Fit functional form to likelihood h(c)
var.flags.fitRates = 1;         % Fit rMin, rMax in h(c)
var.flags.fitPriorFun = 0;      % Fit as yet undetermined function to prior
var.flags.fitG = 0;             % Fit speed-dependence of likelihood g(v)
var.flags.checkFits = 0;        % Check that all the fits returned from fmincon are within
                                    % some margin (var.info.SDThresh) of
                                    % the mean. This might be begging the
                                    % question, so it's left as an option.
var.flags.plotPMFs = 0;         % Plot PMFs for each condition (supp. fig)
var.flags.pmfFitType = 'normal';   % Or 'wbl'. Fit a cumulative of this type to the PMFs
var.flags.normUnits = 'log';    % Or 'linear'. Units over which to normalize integrated prior

%% Define g(v) for each subject, if we're not fitting it as well
if ~var.flags.fitG
    if nStockerSubjects ~= 0
        var.data(1).G = .25;
        var.data(2).G = .2;
    else
        for subject = 1:var.info.nSubjects
            var.data(subject).G = .25;
        end
    end
end

%% Define percent correct for plotting thresholds
var.info.pCor = .82;

%% Organize
var = organizeVar(var);