function var = organizeData(var)
%
% organizeData pulls out the reference and test velocities and contrasts
% for each subject, and gets indices where necessary for use in the
% log-likelihood function
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    % Pull out four variables of interest
    var.data(subject).refV = var.vTrans.transFun(...
                             var.data(subject).rawData(:, var.info.refVCol));
    var.data(subject).refC = var.data(subject).rawData(:, var.info.refCCol);
    var.data(subject).testV = var.vTrans.transFun(...
                              var.data(subject).rawData(:, var.info.testVCol));
    var.data(subject).testC = var.data(subject).rawData(:, var.info.testCCol);
    
    % Collect unique reference speeds and both contrasts
    var.data(subject).refVels = unique(var.data(subject).refV);
    var.data(subject).refContrasts = unique(var.data(subject).refC);
    var.data(subject).testContrasts = unique(var.data(subject).testC);
    
    % Accrue info about uniques and trials
    if subject == 1
        var.info.nTrials = length(var.data(subject).refV);
        var.info.nUniqueRefVels = length(var.data(subject).refVels);
        var.info.nUniqueRefConts = length(var.data(subject).refContrasts);
        var.info.nUniqueContrasts = length(var.data(subject).testContrasts);
    end
    
    % Collect choices, and preallocate resample array
    var.data(subject).T = nan(var.info.nTrials, var.info.nBoots);
    var.data(subject).T(:, 1) = var.data(subject).rawData(:, var.info.choiceCol);
    
    % Get indices of both contrasts for log-likelihood function
    var.data(subject).refCInds = nan(size(var.data(subject).refC));
    var.data(subject).testCInds = var.data(subject).refCInds;
    for contrast = 1:var.info.nUniqueContrasts
        var.data(subject).refCInds(var.data(subject).refC ==...
        var.data(subject).testContrasts(contrast)) = contrast;
    
        var.data(subject).testCInds(var.data(subject).testC ==...
        var.data(subject).testContrasts(contrast)) = contrast;
    end
end

%% Organize
var = organizeVar(var);