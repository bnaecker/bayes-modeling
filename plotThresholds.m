function var = plotThresholds(var)
%
% plotThresholds(var) plots the absolute and relative thresholds for
% subjects performing the speed discrimination in Stocker & Simoncelli 2006
%
% (c) benjamin naecker UT Austin 1 Jun 2011 benjamin.naecker@gmail.com

%% marker colors
threshFC = {'w'; 'k'};
matchFC = {0.*[1 1 1]; .2.*[1 1 1]; .4.*[1 1 1]; .6.*[1 1 1]; .8.*[1 1 1]};

badI = zeros(var.info.nBoots);
%% Go through each subject
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    %% Make one figure for each subject
    var.handles(subject).matchingFig = figure;
    set(var.handles(subject).matchingFig, 'Name', ['Subject ' num2str(subject) ...
        ' Matching Speeds'], 'Color', 'w');
    var.handles(subject).threshFig = figure;
    set(var.handles(subject).threshFig, 'Name', ['Subject ' num2str(subject) ...
        ' Thresholds'], 'Color', 'w');
    
    %% Setup 'var' to receive threshold estimates, one for each strap
    var.stats(subject).matchSpeed = nan(var.info.nUniqueRefConts, ...
        var.info.nUniqueRefVels, var.info.nBoots);
    var.stats(subject).relMatchSpeed = var.stats(subject).matchSpeed;
    var.stats(subject).absThresh = var.stats(subject).matchSpeed;
    var.stats(subject).relThresh = var.stats(subject).matchSpeed;
    
    %% Mean/errbar threshold estimates
    var.stats(subject).meanMatchSpeed = nan(var.info.nUniqueRefConts, ...
        var.info.nUniqueRefVels);
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
            
            % Find the speed that gives the nearest proprotion correct
            % (matching speed, etc) for each individual bootstrap
            for boot = 1:var.info.nBoots
                % From fitted psychometric functions
                if any(isfinite(var.stats(subject).pmf{c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(boot, :)))
                    [~, i] = findnearest(.5, var.stats(subject).pmf{...
                        c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(boot, :));
                    if rem(length(i), 2) == 1
                        i = median(i);
                        badI(boot) = 1;
                    else
                        i = i(length(i)/2);
                        badI(boot) = 1;
                    end
                
                    % Calculate speeds, thresholds
                    var.stats(subject).matchSpeed(c, v, boot) = var.stats(subject).pmfDomain{...
                        c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)}(i);
                    var.stats(subject).relMatchSpeed(c, v, boot) = var.stats(subject).matchSpeed(c, v, boot) / ...
                        var.vTrans.iTransFun(var.data(subject).refVels(v));
                    var.stats(subject).absThresh(c, v, boot) = abs(var.stats(subject).matchSpeed(c, v, boot) - ...
                        var.vTrans.iTransFun(var.data(subject).refVels(v)));
                    var.stats(subject).relThresh(c, v, boot) = var.stats(subject).absThresh(c, v, boot) / ...
                        var.vTrans.iTransFun(var.data(subject).refVels(v));
                else
%                     set(var.handles(subject).pmfAx(c, v, var.data(subject).testContrasts == var.data(subject).refContrasts(c)), 'FontSize', 2);
%                     error('you have a data probrem. subj = %d, refC = %d, refV = %d, testC = %d', subject, c, v, ...
%                         find(var.data(subject).testContrasts == var.data(subject).refContrasts(c)));
%                     var.errList(subject,c,v,var.data(subject).testContrasts == var.data(subject).refContrasts(c)) = true;

                    %%%%%%%%%%%% we're not getting an estimate because
                    %%%%%%%%%%%% fmincon returned NaN's for the params of
                    %%%%%%%%%%%% the cdf fit to the data.
                end
            end
            
            %% Get average matching speeds and thresholds, and SD estimates
            var.stats(subject).meanMatchSpeed(c, v) = mean(var.stats(subject).matchSpeed(c, v, :));
            var.stats(subject).matchSpeedStd(c, v) = std(var.stats(subject).matchSpeed(c, v, :));
            var.stats(subject).meanRelMatchSpeed(c, v) = mean(var.stats(subject).relMatchSpeed(c, v, :));
            var.stats(subject).relMatchSpeedStd(c, v) = std(var.stats(subject).relMatchSpeed(c, v, :));
            var.stats(subject).meanAbsThresh(c, v) = mean(var.stats(subject).absThresh(c, v, :));
            var.stats(subject).absThreshStd(c, v) = std(var.stats(subject).absThresh(c, v, :));
            var.stats(subject).meanRelThresh(c, v) = mean(var.stats(subject).relThresh(c, v, :));
            var.stats(subject).relThreshStd(c, v) = std(var.stats(subject).relThresh(c, v, :));
            
        end
        
        %% First plot relative matching speeds for each contrast
        var.handles(subject).matchingFig;
        var.handles(subject).matchingAx(c) = subplot(1, 2, c); hold on;
        if c == 1
            set(var.handles(subject).matchingAx(c), 'XLim', [.04 15], 'YLim', [.4 2], 'TickDir', 'out', 'XScale', 'log', 'YScale', 'linear');
            var.handles(subject).matchingTitle = title(['Relative Matching Speeds, Ref Contrast = ' num2str(var.data(subject).refContrasts(c))]);
            var.handles(subject).matchingXLab = xlabel('Velocity');
            var.handles(subject).matchingYLab = ylabel('Test V / Ref V');
            var.handles(subject).matchingUnityLine = plot(logspace(log10(.04), log10(15), 1000), ...
                1, 'Color', [.8 .8 .8], 'LineStyle', '-');
        end
        for testC = 1:var.info.nUniqueContrasts
            var.handles(subject).lowCMatchingLine(testC) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
                var.stats(subject).meanRelMatchSpeed(c, :), 'Marker', 'o', 'Color', matchFC{c});
        end
        
        
        %% Plot the absolute thresholds for each reference contrast
        var.handles(subject).threshFig;
        
        var.handles(subject).absThreshAx = subplot(121); hold on;
        if c == 1
            set(var.handles(subject).absThreshAx, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'log', ...
                'XLim', [.04 15]);
            var.handles(subject).absThreshTitle = title('Absolute Thresholds');
            var.handles(subject).absThreshXLab = xlabel('Velocity');
            var.handles(subject).absThreshYLab = ylabel('Absolute Threshold \Delta V');
        end
        var.handles(subject).absThreshLine(c) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
            var.stats(subject).meanAbsThresh(c, :), 'Marker', 'o', ...
            'MarkerFaceColor', threshFC{c}, 'LineStyle', '--', 'Color', 'k');
        
%         % Plot error bars
%         var.handles(subject).absThreshErrs(c, :) = line(...
%             repmat(var.vTrans.iTransFun(var.data(subject).refVels'), 2, 1), ...
%             [var.stats(subject).meanAbsThresh(c, :) + 2 * var.stats(subject).absThreshStd(c, :); ...
%             var.stats(subject).meanAbsThresh(c, :) - 2 * var.stats(subject).absThreshStd(c, :)], ...
%             'LineStyle', '-', 'Color', [.3 .3 .3]);
%         
%         % Turn off legend entries for errorbars
%         for i = 1:length(var.handles(subject).absThreshErrs(c, :))
%             set(get(get(var.handles(subject).absThreshErrs(c, i), 'Annotation'), ...
%                 'LegendInformation'), 'IconDisplayStyle', 'off');
%         end
        
        %% Plot the relative thresholds for each reference contrast
        var.handles(subject).relThreshAx = subplot(122); hold on;
        if c == 1
            set(var.handles(subject).relThreshAx, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'linear', ...
                'XLim', [.04 15]);
            var.handles(subject).relThreshTitle = title('Relative Thresholds');
            var.handles(subject).relThreshXLab = xlabel('Velocity');
            var.handles(subject).relThreshYLab = ylabel('Relative Threshold \Delta V / V');
        end
        var.handles(subject).relThreshLine(c) = plot(var.vTrans.iTransFun(var.data(subject).refVels), ...
            var.stats(subject).meanRelThresh(c, :), 'Marker', 's', ...
            'MarkerFaceColor', threshFC{c}, 'LineStyle', '--', 'Color', 'k');
        
%         % Plot error bars
%         var.handles(subject).relThreshErrs(c, :) = line(...
%             repmat(var.vTrans.iTransFun(var.data(subject).refVels'), 2, 1), ...
%             [var.stats(subject).meanRelThresh(c, :) + 2 * var.stats(subject).relThreshStd(c, :); ...
%             var.stats(subject).meanRelThresh(c, :) - 2 * var.stats(subject).relThreshStd(c, :)], ...
%             'LineStyle', '-', 'Color', [.3 .3 .3]);
%         
%         % Turn off legend entries for errorbars
%         for i = 1:length(var.handles(subject).absThreshErrs(c, :))
%             set(get(get(var.handles(subject).relThreshErrs(c, i), 'Annotation'), ...
%                 'LegendInformation'), 'IconDisplayStyle', 'off');
%         end
        
    end
    % Legends
    var.handles(subject).absThreshLeg = legend(var.handles(subject).absThreshAx, ...
        arrayfun(@num2str, var.data(subject).refContrasts, 'UniformOutput', 0), 'Box', 'off', 'EdgeColor', [1 1 1], 'Location', 'NorthWest');
    set(get(var.handles(subject).absThreshLeg, 'Title'), 'String', 'Contrasts');
    var.handles(subject).relThreshLeg = legend(var.handles(subject).relThreshAx, ...
        arrayfun(@num2str, var.data(subject).refContrasts, 'UniformOutput', 0), 'Box', 'off', 'EdgeColor', [1 1 1], 'Location', 'NorthWest');
    set(get(var.handles(subject).relThreshLeg, 'Title'), 'String', 'Contrasts');
end

%% Organize
var = organizeVar(var);