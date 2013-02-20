function var = plotThresholdsTEST(var)
%
% plotThresholds(var) plots the absolute and relative thresholds for
% subjects performing the speed discrimination in Stocker & Simoncelli 2006
%
% (c) benjamin naecker UT Austin 1 Jun 2011 benjamin.naecker@gmail.com

%% Plot setup
% threshFC = {'w'; 'k'};
matchFC = {1.*[1 1 1]; []; .8.*[1 1 1]; .5.*[1 1 1]; .2.*[1 1 1]; []; 0.*[1 1 1]};
matchXLims = [.04 15];
matchYLims = [.4 2;
             .4 4;
             .01 10;
             0 1];

% badI = zeros(var.info.nBoots);
%% Go through each subject
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    %% Setup 'var' to receive threshold estimates, one for each strap
    var.stats(subject).matchSpeed = cell(var.info.nUniqueRefConts, ...
        var.info.nUniqueRefVels, var.info.nUniqueContrasts);
    var.stats(subject).relMatchSpeed = var.stats(subject).matchSpeed;
    var.stats(subject).absThresh = var.stats(subject).matchSpeed;
    var.stats(subject).relThresh = var.stats(subject).matchSpeed;
    
    %% Mean/errbar threshold estimates
    var.stats(subject).meanMatchSpeed = cell(var.info.nUniqueRefConts, ...
        var.info.nUniqueRefVels, var.info.nUniqueContrasts);
    var.stats(subject).matchSpeedStd = var.stats(subject).meanMatchSpeed;
    var.stats(subject).meanRelMatchSpeed = var.stats(subject).meanMatchSpeed;
    var.stats(subject).relMatchSpeedStd = var.stats(subject).meanMatchSpeed;
    var.stats(subject).meanAbsThresh = var.stats(subject).meanMatchSpeed;
    var.stats(subject).absThreshStd = var.stats(subject).meanMatchSpeed;
    var.stats(subject).meanRelThresh = var.stats(subject).meanMatchSpeed;
    var.stats(subject).relThreshStd = var.stats(subject).meanMatchSpeed;
    
    %% Calculate thresholds
    for c = 1:var.info.nUniqueRefConts
        for v = 1:var.info.nUniqueRefVels
            for testC = 1:var.info.nUniqueContrasts
                    
                % Find speed that gives proportion correct nearest to .5
                % (i.e., matching speed) for each individual bootstrap/PMF
                for boot = 1:var.info.nBoots
                    % as long as the pmf is not empty or all nans here,
                    % calculate
                    if ~isempty(var.stats(subject).pmf{c, v, testC}) && any(isfinite(var.stats(subject).pmf{c, v, testC}(boot, :)))
                        [~, i] = findnearest(.5, var.stats(subject).pmf{c, v, testC}(boot, :));
                        if rem(length(i), 2) == 1 % if there are an odd number of indices that match)
                            i = median(i); % take the median
                        else
                            i = i(length(i)/2); % else take one of the middle ones
                        end
                        
                        % Calculate the speeds and 'thresholds' at this point
                        var.stats(subject).matchSpeed{c, v, testC}(boot) = var.stats(subject).pmfDomain{c, v, testC}(i);
                        var.stats(subject).relMatchSpeed{c, v, testC}(boot) = var.stats(subject).matchSpeed{c, v, testC}(boot) / ...
                            var.vTrans.iTransFun(var.data(subject).refVels(v));
                        var.stats(subject).absThresh{c, v, testC}(boot) = var.stats(subject).matchSpeed{c, v, testC}(boot) - ...
                            var.vTrans.iTransFun(var.data(subject).refVels(v));
                        var.stats(subject).relThresh{c, v, testC}(boot) = var.stats(subject).absThresh{c, v, testC}(boot) / ...
                            var.vTrans.iTransFun(var.data(subject).refVels(v));
                        
                        %  to plot it
                        var.stats(subject).msToPlot{c, v, testC} = var.stats(subject).matchSpeed{c,v,testC}(1);
                        var.stats(subject).rmsToPlot{c, v, testC} = var.stats(subject).relMatchSpeed{c,v,testC}(1);
                    else
                        %%%% Still need to figure out what to do here. This
                        %%%% means that were getting nans in the PMFs,
                        %%%% probably from the fact that some are not being
                        %%%% fit appropriately. WHY?
                    end
                end
                
                % Get the variances of the straps
%                 var.stats(subject).meanMatchSpeed{c, v, testC} = mean(var.stats(subject).matchSpeed{c, v, testC}, 2);
                var.stats(subject).matchSpeedStd{c, v, testC} = std(var.stats(subject).matchSpeed{c, v, testC});
%                 var.stats(subject).meanRelMatchSpeed{c, v, testC} = mean(var.stats(subject).relMatchSpeed{c, v, testC}, 2);
                var.stats(subject).relMatchSpeedStd{c, v, testC} = std(var.stats(subject).relMatchSpeed{c, v, testC});
%                 var.stats(subject).meanAbsThresh{c, v, testC} = mean(var.stats(subject).absThresh{c, v, testC}, 2);
                var.stats(subject).absThreshStd{c, v, testC} = std(var.stats(subject).absThresh{c, v, testC});
%                 var.stats(subject).meanRelThresh{c, v, testC} = mean(var.stats(subject).relThresh{c, v, testC}, 2);
                var.stats(subject).relThreshStd{c, v, testC} = std(var.stats(subject).relThresh{c, v, testC});
            end
        end
    end
    
    %% Plot matching speeds
    
    % Make a figure
    var.handles(subject).matchingFig = figure;
    set(var.handles(subject).matchingFig, 'Name', ['Subject ' num2str(subject) ...
        ' Matching Speeds'], 'Color', 'w');
    
    for c = 1:var.info.nUniqueRefConts
        for testC = 1:var.info.nUniqueContrasts
            
            % if the testC is not one of the refC's, plot it
            if any(testC == find(var.data(subject).testContrasts ~= var.data(subject).refContrasts(1) ...
                & var.data(subject).testContrasts ~= var.data(subject).refContrasts(2)))
            
                % First plot relative matching speeds for each test
                % contrast
                figure(var.handles(subject).matchingFig);
                var.handles(subject).matchingAx(c) = subplot(1, 2, c); hold on;
                set(var.handles(subject).matchingAx(c), 'XLim', matchXLims, 'YLim', matchYLims(c, :), 'TickDir', 'out', 'XScale', 'log', 'YScale', 'linear');
                var.handles(subject).matchingTitle = title(['Relative Matching Speeds, Ref Contrast = ' num2str(var.data(subject).refContrasts(c))]);
                var.handles(subject).matchingXLab = xlabel('Velocity');
                var.handles(subject).matchingYLab = ylabel('Test V / Ref V');
                var.handles(subject).matchingUnityLine = plot(logspace(log10(matchXLims(1)), log10(matchXLims(2)), 1000), ...
                    1, 'Color', [.8 .8 .8], 'LineStyle', '-');
                var.handles(subject).matchLine(c, testC) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
                    [var.stats(subject).rmsToPlot{c, :, testC}], 'Marker', 'o', 'Color', 'k', ...
                    'MarkerFaceColor', matchFC{testC}, 'LineStyle', '--');
            end
        end
    end
    
%     %% Plot thresholds
%     
%     % Make a figure
%     var.handles(subject).threshFig = figure;
%     set(var.handles(subject).threshFig, 'Name', ['Subject ' num2str(subject) ...
%         ' Thresholds'], 'Color', 'w');
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
%         for v = 1:var.info.nUniqueRefVels
%             
%             % Find the speed that gives the nearest proprotion correct
%             % (matching speed, etc) for each individual bootstrap
%             for boot = 1:var.info.nBoots
%                 % From fitted psychometric functions
%                 if any(isfinite(var.stats(subject).pmf{c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(boot, :)))
%                     [~, i] = findnearest(.5, var.stats(subject).pmf{...
%                         c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(boot, :));
%                     if rem(length(i), 2) == 1
%                         i = median(i);
%                         badI(boot) = 1;
%                     else
%                         i = i(length(i)/2);
%                         badI(boot) = 1;
%                     end
%                 
%                     % Calculate speeds, thresholds
%                     var.stats(subject).matchSpeed(c, v, boot) = var.stats(subject).pmfDomain{...
%                         c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(i);
%                     var.stats(subject).relMatchSpeed(c, v, boot) = var.stats(subject).matchSpeed(c, v, boot) / ...
%                         var.vTrans.iTransFun(var.data(subject).refVels(v));
%                     var.stats(subject).absThresh(c, v, boot) = abs(var.stats(subject).matchSpeed(c, v, boot) - ...
%                         var.vTrans.iTransFun(var.data(subject).refVels(v)));
%                     var.stats(subject).relThresh(c, v, boot) = var.stats(subject).absThresh(c, v, boot) / ...
%                         var.vTrans.iTransFun(var.data(subject).refVels(v));
%                 else
% %                     set(var.handles(subject).pmfAx(c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)), 'FontSize', 2);
% %                     error('you have a data probrem. subj = %d, refC = %d, refV = %d, testC = %d', subject, c, v, ...
% %                         find(var.data(subject).testContrasts == var.data(subject).refContrasts(c)));
% %                     var.errList(subject,c,v,var.data(subject).testContrasts == var.data(subject).refContrasts(c)) = true;
% 
%                     %%%%%%%%%%%% we're not getting an estimate because
%                     %%%%%%%%%%%% fmincon returned NaN's for the params of
%                     %%%%%%%%%%%% the cdf fit to the data.
%                 end
%             end
%             
%             %% Get average matching speeds and thresholds, and SD estimates
%             var.stats(subject).meanMatchSpeed(c, v) = mean(var.stats(subject).matchSpeed(c, v, :));
%             var.stats(subject).matchSpeedStd(c, v) = std(var.stats(subject).matchSpeed(c, v, :));
%             var.stats(subject).meanRelMatchSpeed(c, v) = mean(var.stats(subject).relMatchSpeed(c, v, :));
%             var.stats(subject).relMatchSpeedStd(c, v) = std(var.stats(subject).relMatchSpeed(c, v, :));
%             var.stats(subject).meanAbsThresh(c, v) = mean(var.stats(subject).absThresh(c, v, :));
%             var.stats(subject).absThreshStd(c, v) = std(var.stats(subject).absThresh(c, v, :));
%             var.stats(subject).meanRelThresh(c, v) = mean(var.stats(subject).relThresh(c, v, :));
%             var.stats(subject).relThreshStd(c, v) = std(var.stats(subject).relThresh(c, v, :));
%             
%         end
        
%         %% First plot relative matching speeds for each contrast
%         var.handles(subject).matchingFig;
%         var.handles(subject).matchingAx(c) = subplot(1, 2, c); hold on;
%         if c == 1
%             set(var.handles(subject).matchingAx(c), 'XLim', [.04 15], 'YLim', [.4 2], 'TickDir', 'out', 'XScale', 'log', 'YScale', 'linear');
%             var.handles(subject).matchingTitle = title(['Relative Matching Speeds, Ref Contrast = ' num2str(var.data(subject).refContrasts(c))]);
%             var.handles(subject).matchingXLab = xlabel('Velocity');
%             var.handles(subject).matchingYLab = ylabel('Test V / Ref V');
%             var.handles(subject).matchingUnityLine = plot(logspace(log10(.04), log10(15), 1000), ...
%                 1, 'Color', [.8 .8 .8], 'LineStyle', '-');
%         end
%         for testC = 1:var.info.nUniqueContrasts
%             var.handles(subject).lowCMatchingLine(testC) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
%                 var.stats(subject).meanRelMatchSpeed(c, :), 'Marker', 'o', 'Color', matchFC{c});
%         end
%         
%         
%         %% Plot the absolute thresholds for each reference contrast
%         var.handles(subject).threshFig;
%         
%         var.handles(subject).absThreshAx = subplot(121); hold on;
%         if c == 1
%             set(var.handles(subject).absThreshAx, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'log', ...
%                 'XLim', [.04 15]);
%             var.handles(subject).absThreshTitle = title('Absolute Thresholds');
%             var.handles(subject).absThreshXLab = xlabel('Velocity');
%             var.handles(subject).absThreshYLab = ylabel('Absolute Threshold \Delta V');
%         end
%         var.handles(subject).absThreshLine(c) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
%             var.stats(subject).meanAbsThresh(c, :), 'Marker', 'o', ...
%             'MarkerFaceColor', threshFC{c}, 'LineStyle', '--', 'Color', 'k');
%         
% %         % Plot error bars
% %         var.handles(subject).absThreshErrs(c, :) = line(...
% %             repmat(var.vTrans.iTransFun(var.data(subject).refVels'), 2, 1), ...
% %             [var.stats(subject).meanAbsThresh(c, :) + 2 * var.stats(subject).absThreshStd(c, :); ...
% %             var.stats(subject).meanAbsThresh(c, :) - 2 * var.stats(subject).absThreshStd(c, :)], ...
% %             'LineStyle', '-', 'Color', [.3 .3 .3]);
% %         
% %         % Turn off legend entries for errorbars
% %         for i = 1:length(var.handles(subject).absThreshErrs(c, :))
% %             set(get(get(var.handles(subject).absThreshErrs(c, i), 'Annotation'), ...
% %                 'LegendInformation'), 'IconDisplayStyle', 'off');
% %         end
%         
%         %% Plot the relative thresholds for each reference contrast
%         var.handles(subject).relThreshAx = subplot(122); hold on;
%         if c == 1
%             set(var.handles(subject).relThreshAx, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'linear', ...
%                 'XLim', [.04 15]);
%             var.handles(subject).relThreshTitle = title('Relative Thresholds');
%             var.handles(subject).relThreshXLab = xlabel('Velocity');
%             var.handles(subject).relThreshYLab = ylabel('Relative Threshold \Delta V / V');
%         end
%         var.handles(subject).relThreshLine(c) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
%             var.stats(subject).meanRelThresh(c, :), 'Marker', 's', ...
%             'MarkerFaceColor', threshFC{c}, 'LineStyle', '--', 'Color', 'k');
%         
% %         % Plot error bars
% %         var.handles(subject).relThreshErrs(c, :) = line(...
% %             repmat(var.vTrans.iTransFun(var.data(subject).refVels'), 2, 1), ...
% %             [var.stats(subject).meanRelThresh(c, :) + 2 * var.stats(subject).relThreshStd(c, :); ...
% %             var.stats(subject).meanRelThresh(c, :) - 2 * var.stats(subject).relThreshStd(c, :)], ...
% %             'LineStyle', '-', 'Color', [.3 .3 .3]);
% %         
% %         % Turn off legend entries for errorbars
% %         for i = 1:length(var.handles(subject).absThreshErrs(c, :))
% %             set(get(get(var.handles(subject).relThreshErrs(c, i), 'Annotation'), ...
% %                 'LegendInformation'), 'IconDisplayStyle', 'off');
% %         end
    % Legends
%     var.handles(subject).absThreshLeg = legend(var.handles(subject).absThreshAx, ...
%         arrayfun(@num2str, var.data(subject).refContrasts, 'UniformOutput', 0), 'Box', 'off', 'EdgeColor', [1 1 1], 'Location', 'NorthWest');
%     set(get(var.handles(subject).absThreshLeg, 'Title'), 'String', 'Contrasts');
%     var.handles(subject).relThreshLeg = legend(var.handles(subject).relThreshAx, ...
%         arrayfun(@num2str, var.data(subject).refContrasts, 'UniformOutput', 0), 'Box', 'off', 'EdgeColor', [1 1 1], 'Location', 'NorthWest');
%     set(get(var.handles(subject).relThreshLeg, 'Title'), 'String', 'Contrasts');
end

%% Organize
% var = organizeVar(var);