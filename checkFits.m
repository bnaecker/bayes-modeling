function [var, good] = checkFits(var)
%
% Check that the values fit to the resampled data are within
% var.info.SDThresh standard deviations. This is sort of a hack, but it
% throws out anything above this many SDs and resamples/refits the
% parameters of the ideal observer.
%
% (c) benjamin naecker UT Austin 27 May 2011 benjamin.naecker@gmail.com

subject = var.info.currentSubject;

%% Check that all values are within the range
hRowMeans = mean(var.params(subject).hHat(2:end, :), 2);
hRealMean = mean(var.params(subject).hHat(1, :));
hMaxSD = max(std(var.params(subject).hHat));
pRowMeans = mean(var.params(subject).prior(2:end, :), 2);
pRealMean = mean(var.params(subject).prior(1, :));
pMaxSD = max(std(var.params(subject).prior));
checkInds = hRowMeans > hRealMean + var.info.SDThresh .* hMaxSD;
bootNums = find([0; checkInds]);

%% If any values are outside the range, refit just those values
if any(checkInds)
    
    % Notify
    fprintf(['\n\nWARNING: Found %d outlier value(s) greater than %d SDs from '...
        'the mean for original data from subject %d. \n'...
        'Resampling and refitting those values.\n'], sum(checkInds), var.info.SDThresh, subject);
    
    conditions = cartprod(var.data(subject).refVels, var.data(subject).refContrasts, ...
        var.data(subject).testContrasts);
    
    boot = 1;
    while boot <= length(bootNums)
        
        % Resample appropriately
        for condition = 1:size(conditions, 1);
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
                var.data(subject).T(conditionInds, bootNums(boot)) = sample;
            end
        end
        
        % Notify which boot is being fitted
        fprintf('Redoing boot %d (%d of %d)...', bootNums(boot), boot, sum(checkInds));
        
        % Redefine LL function because we've updated fiels of var
        var.info.currentBoot = bootNums(boot);
        var.llOpt(subject).llFun = @(prs)(-logli_var(prs, var));
        
        % Run fmincon
        [var.llOpt(subject).prsHat(:, bootNums(boot)) ...
        var.llOpt(subject).fval(bootNums(boot)) ...
        var.llOpt(subject).eflag(bootNums(boot))] = ...
            fmincon(var.llOpt(subject).llFun, var.llOpt(subject).prs0, [], [], [], [], ...
            var.llOpt(subject).LB, var.llOpt(subject).UB, [], var.llOpt(subject).options);
        
        % If the fitting converged, save appropriately, else redo
        if any(var.llOpt(subject).eflag(bootNums(boot)) == [0 -2])
            fprintf('did not converge! Refitting.\n');
        else
            fprintf('converged @ %.2f\n', var.llOpt(subject).fval(bootNums(boot)));
            
            % Save output appropriately
            var.params(subject).hHat(bootNums(boot), :) = ...
                var.llOpt(subject).prsHat(1:7, bootNums(boot));
            var.params(subject).slopeHat(bootNums(boot), :) = ...
                var.llOpt(subject).prsHat(8:end, bootNums(boot));
            
            % Include g(v), if we're not fitting
            var.params(subject).likeWidth = var.params(subject).hHat(1, :) ./ var.data(subject).G;
            
            % Reconstruct prior from these prior slopes
            yOffset = 20;
            pr = cumtrapz(var.data(subject).refVels, -var.params(subject).slopeHat(bootNums(boot), :), 2) + yOffset;
            pr = pr./sum(pr);
            var.params(subject).prior(bootNums(boot), :) = pr;
            
            % Increment boot counter
            boot = boot + 1;
        end
    end
    
    % Recheck that all values are good
    rowMeans = mean(var.params(subject).hHat(2:end, :), 2);
    realMean = mean(var.params(subject).hHat(1, :));
    maxSD = max(std(var.params(subject).hHat));
    checkInds = rowMeans > realMean + var.info.SDThresh .* maxSD;
    if sum(checkInds) == 0
        good = 1;
    else
        good = 0;
    end
else
    good = 1;
end

%% Organize
var = organizeVar(var);