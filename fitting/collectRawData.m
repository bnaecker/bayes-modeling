function ds = collectRawData(ds)
%
% FUNCTION ds = collectRawData(ds)
%
% This function pulls out psychophysical discrimination information
% relevant for fitting the ideal observer model from the raw data files.
% This includes reference and test velocities and contrasts, as well as
% indices to these values for logical indexing during the fitting
% procedure.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

for si = 1:ds.info.nSubjects
    fprintf('Collecting velocities, contrasts, and choices for subject %d ... ', si);
    ds.info.currentSubject = si;
    
    % pull out contrast and velocity for reference and test stimuli
    ds.data(si).refV = ds.data(si).rawData(:, ds.info.refVCol);
    ds.data(si).refC = ds.data(si).rawData(:, ds.info.refCCol);
    ds.data(si).testV = ds.data(si).rawData(:, ds.info.testVCol);
    ds.data(si).testC = ds.data(si).rawData(:, ds.info.testCCol);
    
    % NOTE: In the data collected by Stocker & Simoncelli, there is a
    % single negative velocity value. They used an adaptive staircase
    % procedure to determine the velocity on the next trial, and such a
    % procedure can of course produce a negative value. Though why this was
    % not censored (as the gratings would then be drifting in opposite
    % directions and discrimination should be relatively easy) is unclear.
    % The authors make no mention of this in their paper or methods
    % section. However, problems arise when trying to fit a cumulative
    % Weibull psychometric function. In this case, the function will return
    % an imaginary value, in which case most fitting toolboxes (including
    % MATLAB's) will throw an error. How this was remedied in the authors'
    % analysis is not specified, but here we simply take the absolute value
    % of the velocity.
    ds.data(si).testV = abs(ds.data(si).testV);
    
    % log-transform real data (not simulated);
    if ~ds.info.isSimData
        ds.data(si).refV = ds.velTrans.transFun(ds.data(si).refV);
        ds.data(si).testV = ds.velTrans.transFun(ds.data(si).testV);
    end
    
    % collect unique reference speeds and both contrasts
    ds.data(si).refVels = unique(ds.data(si).refV);
    ds.data(si).refContrasts = unique(ds.data(si).refC);
    ds.data(si).testContrasts = unique(ds.data(si).testC);
    
    % accrue information about uniques and trials
    if si == 1
        ds.info.nTrials = length(ds.data(si).refV);
        ds.info.nUniqueRefVels = length(ds.data(si).refVels);
        ds.info.nUniqueRefContrasts = length(ds.data(si).refContrasts);
        ds.info.nUniqueContrasts = length(ds.data(si).testContrasts);
    end
    
    % collect choices and preallocate resample arrays
    ds.data(si).T = nan(ds.info.nTrials, ds.info.nBoots);
    ds.data(si).T(:, 1) = ds.data(si).rawData(:, ds.info.choiceCol);
    
    % get indices of both contrasts for arrangement during minimization of
    % the log-likelihood of the model
    ds.data(si).refCInds = nan(size(ds.data(si).refC));
    ds.data(si).testCInds = ds.data(si).refCInds;
    for ci = 1:ds.info.nUniqueContrasts
        ds.data(si).refCInds(ds.data(si).refC == ...
            ds.data(si).testContrasts(ci)) = ci;
        ds.data(si).testCInds(ds.data(si).testC == ...
            ds.data(si).testContrasts(ci)) = ci;
    end
    fprintf('done.\n');
end

%% organize
ds = orderDataStruct(ds);