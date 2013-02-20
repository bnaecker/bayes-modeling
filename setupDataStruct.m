function ds = setupDataStruct(args)
%
% FUNCTION ds = setupDataStruct(args)
%
% Loads the requested files from disk, which contain psychophysical
% discrimination data. The function assumes that the data is contained in a
% folder within the current working directory, which is called 'Data', and
% that the data is in the format described in the README.
%
% INPUT:    args - the parameter-values pairs output of
%               parseSpeedDiscriminationInput.m
%
% OUTPUT:   ds - the main data structure, which from here on is the only
%               input and output to all functions.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

%% intialize the data structure
ds = struct(...
    'data', {[]}, ...           % data from each subject
    'flags', {[]}, ...          % flags for fitting, plotting, etc.
    'likeOpt', {[]}, ...        % optimization setup for fitting likelihood
    'handles', {[]}, ...        % all graphics handles
    'info', {[]}, ...           % general information about the model/data
    'llOpt', {[]}, ...          % optimization setup for log-likelihood of the model
    'params', {[]}, ...         % actual parameters returned by the optimization
    'stats', {[]}, ...          % statistics of the resampled fits
    'pmfs', {[]}, ...           % pscyhometric functions, model-predicted and fitted
    'velTrans', {[]} ...        % info for the log-transform of the velocity domain
    );

%% setup the file structure and load the data
fprintf('\nFinding files and loading data ... ');
ds.info.baseDir = fileparts(pwd);
ds.info.dataDir = fullfile(ds.info.baseDir, 'Data');
if args.nSubjects == 0
    ds.info.nSubjects = 1;
    if strcmp(args.priorType, 'gaussian');
        d = load(fullfile(ds.info.dataDir, 'simData_gaussian.mat'));
        field = 'm';
    elseif strcmp(args.priorType, 'loglinear');
        d = load(fullfile(ds.info.dataDir, 'simData_loglinear.mat'));
        field = 'm';
    else
        error('runSpeedDiscriminationModel:setupDataStruct:unknownPriorType', ...
            'Supported prior types are "loglinear" and "gaussian"');
    end
    ds.data(1).rawData = d.(field);
%     ds.data.sigNse = [d.sig1 d.sig2];
    ds.info.isSimData = true;
else
    ds.info.nSubjects = args.nSubjects;
    ds.info.isSimData = false;
    ds.info.dataFiles = cell(ds.info.nSubjects, 1);
    for si = 1:ds.info.nSubjects
        ds.info.dataFiles{si} = ['s' num2str(si) '.mat'];
        d = load(fullfile(ds.info.dataDir, ds.info.dataFiles{si}));
        fieldName = ['s' num2str(si)];
        ds.data(si).rawData = d.(fieldName)';
    end
end
fprintf('done.\nSetting up parameters and flags ... ');

%% definition of column indices
ds.info.refVCol = 1;        % speed of reference stimulus
ds.info.refCCol = 2;        % contrast of reference stimulus
ds.info.testVCol = 3;       % speed of test stimulus
ds.info.testCCol = 4;       % contrast of test stimulus
ds.info.stairCol = 5;       % indicates which of two interleaved staircases
ds.info.dirCol = 6;         % direction of stimuli
ds.info.placeCol = 7;       % location of reference stimulus (1 == right)
ds.info.rtCol = 8;          % reaction time
ds.info.choiceCol = 9;      % subject's choice (1 == test seen faster)
ds.info.adaptCol = 10;      % step direction of staircase
ds.info.correctCol = 11;    % correct choice

%% setup log-transformation of velocity domain
ds.velTrans.v0 = 0.3;
ds.velTrans.transFun = ...
    @(v) (log(1 + (v ./ ds.velTrans.v0)));      % transformation to log-axis
ds.velTrans.iTransFun = ...
    @(v) (expm1(v) .* ds.velTrans.v0);          % inverse log-transformation

%% flags
ds.flags = args;
ds.flags = rmfield(ds.flags, {'nSubjects', 'nBoots', 'nSEs'});
ds.info.nSEs = args.nSEs;
% if ds.info.isSimData
%     ds.info.nBoots = 1;
% else
    ds.info.nBoots = args.nBoots;
% end

%% define speed-dependence of likelihood widths
if ds.flags.fitLikeSpeed
    % placeholder to include functionality if desired.
else
    if ds.info.isSimData
        ds.data.G = 1;
    else
        % these values are simply read off figure 4 from Stocker & Simoncelli
        % 2006
        ds.data(1).G = 0.25;
        if ds.info.nSubjects > 1
            ds.data(2).G = 0.2;
        end
    end
end

%% percent correct for determining tresholds
ds.info.pCor = 0.82;

%% organize
ds = orderDataStruct(ds);
fprintf('done.\n');