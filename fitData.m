function var = fitData(var)
%
% fitData returns the prior and likelihood h(c) that maximizes the
% log-likelihood function
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

% Notify
fprintf('\nRunning optimization to fit prior and likelihood widths.\n');
[var.llOpt.numFails] = deal(zeros(var.info.nBoots, 1));

% Run on each subject
for subject = 1:var.info.nSubjects
    
    % Notify caller the subject being analyzed
    fprintf('\nSubject %d\n', subject);
    
    % Setup indices
    var.info.currentSubject = subject;
    var.info.currentBoot = 1;
    
    % Define all 3-conditions to index into choices for resampling
    conditions = cartprod(var.data(subject).refVels, var.data(subject).refContrasts, ...
                          var.data(subject).testContrasts);
                      
    while var.info.currentBoot <= var.info.nBoots
        
        % Resample appropriately
        if var.info.currentBoot > 1 
            for condition = 1:size(conditions, 1)
                currentCondition = conditions(condition, :);
                conditionInds = var.data(subject).refV == currentCondition(1) & ...
                                var.data(subject).refC == currentCondition(2) & ...
                                var.data(subject).testC == currentCondition(3);
                if sum(conditionInds) > 0
                    uniqueV = unique(var.data(subject).testV(conditionInds));
                    dat = var.data(subject).T(conditionInds, 1);
                    sample = nan(size(dat));
                    for v = 1:length(uniqueV)
                        vInds = uniqueV(v) == var.data(subject).testV(conditionInds);
                        if sum(vInds) == 1 
                            sample(vInds) = dat(vInds);
                        else 
                            sample(vInds) = randsample(dat(vInds), sum(vInds), true)';
                        end
                    end
                    var.data(subject).T(conditionInds, var.info.currentBoot) = sample;
                end
            end
        end
        
        % Notify which sample being fitted
        fprintf('Boot %d of %d...', var.info.currentBoot, var.info.nBoots);
        
        % Redefine LL function because we've updated fields of var (boot
        % counter)
        var.llOpt(subject).llFun = @(prs)(-logli_var(prs, var));
        
        % Run fmincon
        [var.llOpt(subject).prsHat(:, var.info.currentBoot) ...
         var.llOpt(subject).fval(var.info.currentBoot) ...
         var.llOpt(subject).eflag(var.info.currentBoot)] = ...
            fmincon(var.llOpt(subject).llFun, var.llOpt(subject).prs0, [], [], [], [], ...
            var.llOpt(subject).LB, var.llOpt(subject).UB, [], var.llOpt(subject).options);
                
        % If the fitting converged, save appropriately, else redo
        switch var.llOpt(subject).eflag(var.info.currentBoot)
            case 0
                var.llOpt.numFails(var.info.currentBoot) = ...
                    var.llOpt.numFails(var.info.currentBoot) + 1;
                if var.llOpt.numFails(var.info.currentBoot) < var.llOpt.numAllowedFails
                    fprintf(['did not converge! MaxFunEvals (%d) or MaxIter (%d) reached, ' ...
                        'consider increasing. Refitting.\n'], var.llOpt(subject).options.MaxFunEvals, ...
                        var.llOpt(subject).options.MaxIter);
                else
                    fprintf('did not converge! \nFitting has now failed %d times!', ...
                        var.llOpt.numFails(var.info.currentBoot));
                    procede = input(' Would you like to procede? [0] or [1]');
                    if ~procede
                        var.error = sprintf(['\nYou elected to quit because' ...
                            'fmincon does not like you. You should investigate' ...
                            'why it is taking more than %d to reach a local minimum.'], ...
                            var.llOpt.options.MaxFunEvals);
                        break
                    end
                end
            case -2
                var.llOpt.numFails(var.info.currentBoot) = ...
                    var.llOpt.numFails(var.info.currentBoot) + 1;
                if numFails >= 10
                else
                    fprintf(['did not converge! No feasible point, ' ...
                        'consider including feasible points. Refitting.\n']);
                end
            otherwise
                fprintf('converged @ %.2f\n', var.llOpt(subject).fval(var.info.currentBoot));
                
                % Save output appropriately
                var.params(subject).hHat(var.info.currentBoot, :) = ...
                    var.llOpt(subject).prsHat(1:var.info.nUniqueContrasts, var.info.currentBoot);
                var.params(subject).slopeHat(var.info.currentBoot, :) = ...
                    var.llOpt(subject).prsHat(var.info.nUniqueContrasts+1:end, var.info.currentBoot);
                
                % Include g(v), if we're not fitting
                var.params(subject).likeWidth = var.params(subject).hHat(1, :) ./ var.data(subject).G;
                
                % Reconstruct prior from these prior slopes
                var = reconstructPrior(var);
%                 yOffset = 20;
%                 pr = cumtrapz(var.data(subject).refVels, -var.params(subject).slopeHat(var.info.currentBoot, :), 2) + yOffset;
%                 pr = pr./sum(pr);
%                 var.params(subject).prior(var.info.currentBoot, :) = pr;
                
                % Increment boot counter
                var.info.currentBoot = var.info.currentBoot + 1;
        end
    end
end

% Notify
if exist('procede', 'var')
    if procede == 0
        fprintf('\nOptimization aborted.\n');
    end
else
    fprintf('\nOptimization completed.\n\n');
end

%% Organize
var = organizeVar(var);