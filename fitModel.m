function ds = fitModel(ds)
%
% FUNCTION ds = fitModel(ds)
%
% The function fitModel.m fits the model assuming the requested prior type
% to the psychophysical data in the data structure ds.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

%% Notify
fprintf('\nRunning optimization to fit prior and likelihood.\n')

%% Run on each subject
for subject = 1:ds.info.nSubjects
    
    %% Notify for each subject
    fprintf('\nSubject %d of %d\n', subject, ds.info.nSubjects);
    
    %% Setup indices
    ds.info.currentSubject = subject;
    ds.info.currentBoot = 1;
    
    %% Define all triplet conditions (refV, refC, testC) to index into
    % choices for resampling
    conditions = cartprod(ds.data(subject).refVels, ...
        ds.data(subject).refContrasts, ...
        ds.data(subject).testContrasts);
    
    % Outer fitting loop (loop over bootstraps)
    while ds.info.currentBoot <= ds.info.nBoots
        
        %% Resampling
        % Check that we're resampling (or just taking the true choices)
        if ds.info.currentBoot > 1
            for ci = 1:size(conditions, 1)
                % Get the current condition triplet and its index
                currentCondition = conditions(ci, :);
                conditionIndices = ...
                    ds.data(subject).refV == currentCondition(1) & ...
                    ds.data(subject).refC == currentCondition(2) & ...
                    ds.data(subject).testC == currentCondition(3);
                
                % Check that there are any trials of this triplet
                if sum(conditionIndices) > 0
                    % Find the unique velocities from this triplet
                    uniqueVels = unique(ds.data(subject).testV(conditionIndices));
                    % Get the subject's choices from this triplet
                    d = ds.data(subject).T(conditionIndices, 1);
                    % Preallocate resample array
                    reSample = NaN(size(d));
                    % For each unique velocity in this triplet
                    for vi = 1:length(uniqueVels)
                        % Find the indices of these velocities
                        vIndices = uniqueVels(vi) == ...
                            ds.data(subject).testV(conditionIndices);
                        % If there is only one
                        if sum(vIndices) == 1
                            % Make the resample that one
                            reSample(vIndices) = d(vIndices);
                        else
                            % Otherwise, randomly resample the data with
                            % replacement
                            reSample(vIndices) = ...
                                randsample(d(vIndices), sum(vIndices), true)';
                        end
                    end
                    % Allocate the resampled data to the choice array T
                    ds.data(subject).T(...
                        conditionIndices, ds.info.currentBoot) = reSample;
                end
            end
        end
        
        %% Notify of which sample being fitted
        fprintf('Boot %d of %d...', ds.info.currentBoot, ds.info.nBoots);
        
        %% Redefine handle to log-likelihood function
        % this must be done since we've updated the data stucture 'ds',
        % which is an argument to the function itself.
        if strcmp(ds.flags.priorType, 'loglinear')
            ds.llOpt(subject).llFun = ...
                @(prs) (-logli_loglinear(prs, ds));
        else
            ds.llOpt(subject).llFun = ...
                @(prs) (-logli_gaussian(prs, ds));
        end
        
        %% Run optimization
        [ds.llOpt(subject).prsHat(:, ds.info.currentBoot) ...
         ds.llOpt(subject).fval(ds.info.currentBoot) ...
         ds.llOpt(subject).eflag(ds.info.currentBoot)] = ...
            fmincon(ds.llOpt(subject).llFun, ds.llOpt(subject).prs0, ...
            [], [], [], [], ds.llOpt(subject).LB, ds.llOpt(subject).UB, ...
            [], ds.llOpt(subject).options);
        
        %% If the fitting converged, save appropriately
        % Note that the data is not saved in the case that the fitting does
        % not converge. Each users must decide what is to be done in these
        % cases.
        switch ds.llOpt(subject).eflag(ds.info.currentBoot)
            case 0      % Number of iterations exceeded MaxIter or MaxFunEvals
            case -2     % No feasible point found
            otherwise   % Converged
                % Notify
                fprintf('converged @ %.2f\n', ds.llOpt(subject).fval(ds.info.currentBoot));
                
                % Save appropriately
                ds.params(subject).hHat(ds.info.currentBoot, :) = ...
                    ds.llOpt(subject).prsHat(1:ds.info.nUniqueContrasts, ...
                    ds.info.currentBoot);
                if strcmp(ds.flags.priorType, 'loglinear')
                    ds.params(subject).slopeHat(ds.info.currentBoot, :) = ...
                        ds.llOpt(subject).prsHat(ds.info.nUniqueContrasts + 1 : end, ...
                        ds.info.currentBoot);
                else
                    ds.params(subject).gamma(ds.info.currentBoot) = ...
                        ds.llOpt(subject).prsHat(end, ds.info.currentBoot); 
                end
                
                % Include g(v) if it is not being fit
                if ds.flags.fitG
                else
                    ds.params(subject).likeWidth = ds.params(subject).hHat ./ ...
                        ds.data(subject).G;
                end
                
                %% Reconstruct the prior from the fitted prior slopes
                ds = reconstructPrior(ds);
                
                %% Increment boot counter
                ds.info.currentBoot = ds.info.currentBoot + 1;
        end
    end
end

%% Organize
ds = orderDataStruct(ds);