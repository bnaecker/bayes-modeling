function ds = collectData(ds)
%
% FUNCTION ds = collectData(ds)
%
% INPUT:    ds = the structure output of setupData
%
% OUTPUT:   ds = same
%
% The function collectData.m pulls out the psychophysical information most
% relevant for fitting the model, such as the various reference and test
% velocities and contrasts, as well as indices to these values for fast
% logical indexing during the fitting procedure.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

fprintf('\nCollecting velocities, contrasts, and choices for ');
for subject = 1:ds.info.nSubjects
    fprintf('subject %d...', subject);
    ds.info.currentSubject = subject;
    
    % Pull out the contrast and velocity for the reference and test stimuli
    ds.data(subject).refV = ds.data(subject).rawData(:, ds.info.refVCol);
    ds.data(subject).refC = ds.data(subject).rawData(:, ds.info.refCCol);
    ds.data(subject).testV = ds.data(subject).rawData(:, ds.info.testVCol);
    ds.data(subject).testC = ds.data(subject).rawData(:, ds.info.testCCol);
    
    % Transform if real data
    if ~ds.info.isSimData
        ds.data(subject).refV = ds.vTrans.transFun(ds.data(subject).refV);
        ds.data(subject).testV = ds.vTrans.transFun(ds.data(subject).testV);
    end
    
    % Collect the unique reference speeds and both contrasts
    ds.data(subject).refVels = unique(ds.data(subject).refV);
    ds.data(subject).refContrasts = unique(ds.data(subject).refC);
    ds.data(subject).testContrasts = unique(ds.data(subject).testC);
    
    % Accrue information about uniques and trials
    if subject == 1
        ds.info.nTrials = length(ds.data(subject).refV);
        ds.info.nUniqueRefVels = length(ds.data(subject).refVels);
        ds.info.nUniqueRefContrasts = length(ds.data(subject).refContrasts);
        ds.info.nUniqueContrasts = length(ds.data(subject).testContrasts);
    end
    
    % Collect choices and preallocate resample arrays
    ds.data(subject).T = NaN(ds.info.nTrials, ds.info.nBoots);
    ds.data(subject).T(:, 1) = ds.data(subject).rawData(:, ds.info.choiceCol);
    
    % Get indices of both contrasts for arrangement during log-likelihood
    % optimization
    ds.data(subject).refCInds = NaN(size(ds.data(subject).refC));
    ds.data(subject).testCInds = ds.data(subject).refCInds;
    for contrast = 1:ds.info.nUniqueContrasts
        ds.data(subject).refCInds(ds.data(subject).refC == ...
            ds.data(subject).testContrasts(contrast)) = contrast;
        ds.data(subject).testCInds(ds.data(subject).testC == ...
            ds.data(subject).testContrasts(contrast)) = contrast;
    end
    fprintf('done. ');
end

%% Organize
ds = orderDataStruct(ds);