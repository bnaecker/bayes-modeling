function var = plotPMFsTEST(var)
%
% Plot the psychometric functions for each 3-condition in the speed
% discrimination task. Cf. online supplementary figure.
%
% (c) benjamin naecker UT Austin 31 May 2011 benjamin.naecker@gmail.com

%% Go through each subject and calculate each PMF
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    %% Plot PMFs if requested
    % Setup cell array to receive variable proportions
    var.stats(subject).propTestFaster = cell(var.info.nUniqueRefConts, ...
        var.info.nUniqueRefVels, var.info.nUniqueContrasts);
    var.stats(subject).meanProps = var.stats(subject).propTestFaster;
    var.stats(subject).propSDs = var.stats(subject).meanProps;
    
    % Loop through each reference contrast
    for refContrastNum = 1:var.info.nUniqueRefConts
        if var.flags.plotPMFs
            % Create figure for each subject and reference contrast
            var.handles(subject).pmfFig= figure;
            set(gcf, 'Name', ['Subject ' num2str(subject) ' PMFs, Ref Contrast = '...
                num2str(var.data(subject).refContrasts(refContrastNum))], 'Color', 'w',...
                'Units', 'Normalized');
            
            % Create labels for direction of increasing ref speed and test
            % contrast
            var.handles(subject).pmfTextAx = axes('Position', [.05 .95 .05 .05], ...
                'Units', 'Normalized', 'Visible', 'off'); hold on;
            var.handles(subject).pmfRefVelTxt = text(0, .2, 'Ref Vel.',...
                'FontSize', 12, 'FontName', 'Helvetica');
            var.handles(subject).pmfRefVelArrow = text(.5, 0, '\downarrow', ...
                'FontSize', 12', 'FontName', 'Helvetica');
            var.handles(subject).pmfTestCArrow = text(0, .6, 'Test Contr. \rightarrow',...
                'FontSize', 12, 'FontName', 'Helvetica');
        end
        
        % Set index to subplots
        subIndex = 0;
        
        % Go through each condition (ref vel & test contr), and plot the
        % number of times the subject indicated the "test" grating was
        % seen to move faster as a function of the test grating speed (i.e., a
        % psychometric function)
        for refVelocityNum = 1:var.info.nUniqueRefVels
            for testContrastNum = 1:var.info.nUniqueContrasts
                
                % Indices to this condition
                refCInds = var.data(subject).refC == var.data(subject).refContrasts(refContrastNum);
                refVInds = var.data(subject).refV == var.data(subject).refVels(refVelocityNum);
                testCInds = var.data(subject).testC == var.data(subject).testContrasts(testContrastNum);
                
                % Get test velocities at this condition
                testVels = unique(var.data(subject).testV(refCInds & refVInds & testCInds,:));
                
                % Check that we have data
                if ~isempty(testVels)
                    
                    % Bootstrap an entire range of the test velocities, so
                    % that we can fit a PMF to it at the same time
                    for boot = 1:var.info.nBoots
                        for testVelocityNum = 1:length(unique(testVels))
                            
                            % Indices again
                            INDS = refCInds & refVInds & testCInds & var.data(subject).testV == ...
                                testVels(testVelocityNum);
                            
                            % Calculate P(testFaster) for each testVel
                            % choie == 1 means reference was faster, so
                            % find choice == 0 to get test faster
                            var.stats(subject).propTestFaster{refContrastNum, refVelocityNum, ...
                                testContrastNum}(boot, testVelocityNum) = sum(var.data(subject).T(INDS, boot)) ./ ...
                                length(var.data(subject).T(INDS, boot));
                        end
                        
                        % Fit desired cumulative distribution to each
                        % bootstrapped set of the entire range of testVels
                        prs0 = [1; 3];
%                         prs0 = [var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum)); ...
%                             1/var.data(subject).testContrasts(testContrastNum)];
                        LB = [1e-15; 1e-15];
                        UB = [100; 100];
                        fitFun = @(prs)(fitCDF(prs, var.flags.pmfFitType, ...
                            var.vTrans.iTransFun(testVels'), var.stats(subject).propTestFaster{...
                            refContrastNum, refVelocityNum, testContrastNum}(boot, :)));
                        var.stats(subject).pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, :) = ...
                            fmincon(fitFun, prs0, [], [], [], [], LB, UB, [], var.hOpt(subject).options);
                        var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, testContrastNum} = ...
                            linspace(.5 * var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum)), ...
                            2* var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum)), 1000);
                        var.stats(subject).pmf{refContrastNum, refVelocityNum, testContrastNum}(boot, :) = ...
                            cdf(var.flags.pmfFitType, var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}, ...
                            var.stats(subject).pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, 1), ...
                            var.stats(subject).pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, 2));
                    end
                    
                    % Calculate mean/SDs of all the bootstrapped choice sets
                    var.stats(subject).meanProps{refContrastNum, refVelocityNum, testContrastNum} = ...
                        mean(var.stats(subject).propTestFaster{refContrastNum, refVelocityNum, ...
                        testContrastNum}, 1);
                    var.stats(subject).propSDs{refContrastNum, refVelocityNum, testContrastNum} = ...
                        std(var.stats(subject).propTestFaster{refContrastNum, refVelocityNum, ...
                        testContrastNum}, 1);
                    
                    % Calculate mean/SDs of the bootstrapped PMF fits
                    var.stats(subject).meanPmf{refContrastNum, refVelocityNum, testContrastNum} = ...
                        mean(var.stats(subject).pmf{refContrastNum, refVelocityNum, testContrastNum}, 1);
                    
%                     % Fit cumulative distribution to the data (transformed
%                     % back to linear domain)
%                     prs0 = [var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum));...
%                         1/var.data(subject).testContrasts(testContrastNum)];
%                     LB = [1e-15; 1e-15];
%                     UB = [100; 100];
%                     fitFun = @(prs)(fitCDF(prs, var.flags.pmfFitType, ...
%                         var.vTrans.iTransFun(testVels'), var.stats(subject).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}(1, :)));
%                     var.stats(subject).pmfPrs(refContrastNum, refVelocityNum, testContrastNum, :) = ...
%                         fmincon(fitFun, prs0, [], [], [], [], LB, UB, [], var.hOpt(subject).options);
%                     var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, testContrastNum} = ...
%                         linspace(.5 * var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum)),...
%                         2 * var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum)), 1000);
%                     var.stats(subject).pmf{refContrastNum, refVelocityNum, testContrastNum} = ...
%                         cdf(var.flags.pmfFitType, var.stats(subject).pmfDomain{refContrastNum, ...
%                         refVelocityNum, testContrastNum}, var.stats(subject).pmfPrs(refContrastNum, ...
%                         refVelocityNum, testContrastNum, 1), var.stats(subject).pmfPrs(refContrastNum, ...
%                         refVelocityNum, testContrastNum, 2));
                end
                
                if var.flags.plotPMFs
                    % Plot into correct subplot
                    subIndex = subIndex + 1;
                    var.handles(subject).pmfAx(refContrastNum, refVelocityNum, testContrastNum) = ...
                        subplot(var.info.nUniqueRefVels, var.info.nUniqueContrasts, subIndex); hold on;
                    
                    % Label the rows/cols
                    if mod(subIndex, var.info.nUniqueContrasts) == 1
                        var.handles(subject).pmfVLabels(refContrastNum, refVelocityNum) = ...
                            ylabel([num2str(var.vTrans.iTransFun(var.data(subject).refVels(refVelocityNum))) ' [deg s^{-1}]      '], ...
                            'FontSize', 12, 'FontWeight', 'Bold', 'Rotation', 0);
                        if subIndex == (var.info.nUniqueRefVels - 1) * var.info.nUniqueContrasts + 1
                            var.handles(subject).pmfCLabels(refContrastNum, testContrastNum) = xlabel(...
                                num2str(var.data(subject).testContrasts(testContrastNum)), ...
                                'FontSize', 12, 'FontWeight', 'Bold');
                        end
                    elseif any(subIndex == (var.info.nUniqueRefVels - 1) * var.info.nUniqueContrasts + 1 : ...
                            var.info.nUniqueRefVels * var.info.nUniqueContrasts)
                        var.handles(subject).pmfVLabels(refContrastNum, testContrastNum) = xlabel(...
                            num2str(var.data(subject).testContrasts(testContrastNum)), ...
                            'FontSize', 12, 'FontWeight', 'Bold');
                    end
                    if ~isempty(testVels)
                        set(var.handles(subject).pmfAx(refContrastNum, refVelocityNum, testContrastNum), ...
                            'TickDir', 'out', 'XLim', [var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, ...
                            testContrastNum}(1) var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}(end)], ...
                            'Units', 'Normalized');
                    end
                    
                    % Plot data and fit
                    if ~isempty(testVels)
                        var.handles(subject).pmfData(refContrastNum, refVelocityNum, testContrastNum) = ...
                            plot(var.vTrans.iTransFun(testVels'), var.stats(subject).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}(1, :), ...
                            'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 5);
                        var.handles(subject).pmfFit(refContrastNum, refVelocityNum, testContrastNum) = ...
                            plot(var.stats(subject).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}, var.stats(subject).pmf{...
                            refContrastNum, refVelocityNum, testContrastNum}(1, :), ...
                            'LineStyle', '-.', 'LineWidth', 2, 'Marker', 'none', 'Color', [1 .3 .3]);
                    else
                        set(var.handles(subject).pmfAx(refContrastNum, refVelocityNum, testContrastNum), ...
                            'Visible', 'off');
                        var.handles(subject).emptyText{refContrastNum, refVelocityNum, testContrastNum} = ...
                            text(.5, .5, 'No Data');
                    end
                end
            end
        end
    end
end

%% Organize
var = organizeVar(var);