function var = plotFits(var)
%
% Plots the likelihood widths and prior for each subject.
%
% (c) benjamin naecker UT Austin 27 May 2011 benjamin.naecker@gmail.com

%% Go through each subject
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    %% Make a figure and a subplot for likelihood and prior
    var.handles(subject).fig = figure;
    set(var.handles(subject).fig, 'Name', ['Subject ' num2str(subject)], ...
        'Color', 'w');
    var.handles(subject).likeAx = subplot(121); hold on; axis square;
    var.handles(subject).likeTitle = title('Likelihood h(c)');
    var.handles(subject).likeXLab = xlabel('Contrast');
    var.handles(subject).priorAx = subplot(122); hold on; axis square;
    var.handles(subject).priorTitle = title('Prior');
    var.handles(subject).priorXLab = xlabel('Speed (deg s^{-1})');
    
    %% Pretty up the axes
    set(var.handles(subject).likeAx, 'XLim', [.04 1], 'YLim', [.5 3], 'XScale', 'log', ...
        'TickDir', 'out', 'FontSize', 12, 'FontName', 'Helvetica');
    set(var.handles(subject).priorAx, 'XLim', [.4 13], 'YLim', [1e-14 2], 'XScale', 'log', ...
        'YScale', 'log', 'TickDir', 'out', 'FontSize', 12, 'FontName', 'Helvetica');
    
    %% Plot likelihood and, if requested, SEs and h(c)
    set(var.handles(subject).fig, 'CurrentAxes', var.handles(subject).likeAx);
    var.handles(subject).likeLine = ...
        plot(var.data(subject).testContrasts, var.params(subject).likeWidth, ...
        'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    if var.flags.plotSmears
        var.handles(subject).likeSEs = ...
            fill([var.data(subject).testContrasts; flipud(var.data(subject).testContrasts)], ...
            [var.params(subject).likeLB'; flipud(var.params(subject).likeUB')], ...
            [.75 .75 .75], 'FaceAlpha', .5, 'EdgeColor', 'none');
    end
    if var.flags.fitLikeFun
        var.handles(subject).hLine = ...
            plot(var.hOpt(subject).domain, var.data(subject).hFunVals, ...
            'LineStyle', '--', 'Marker', 'none', 'Color', 'k');
    end
    
    %% Plot priors and, if requested, SEs
    set(var.handles(subject).fig, 'CurrentAxes', var.handles(subject).priorAx);
    var.handles(subject).priorLine = ...
        plot(var.vTrans.iTransFun(var.data(subject).refVels), var.params(subject).prior(1, :), ...
        'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    if var.flags.plotSmears
        var.handles(subject).priorSEs = ...
            fill(var.vTrans.iTransFun([var.data(subject).refVels; flipud(var.data(subject).refVels)]), ...
            [var.params(subject).prior(1, :)'; flipud(var.params(subject).prUB(1, :)')], ...
            [.75 .75 .75], 'FaceAlpha', .5, 'EdgeColor', 'none');
    end
    
    %% Plot all resamples, if requested
    if var.flags.plotStraps
        var.handles(subject).extraFig = figure;
        set(var.handles(subject).extraFig, 'Name', ['Subject ' num2str(subject) ', Individual Resamples'], ...
            'Color', 'w');
        var.handles(subject).extraLikeAx = subplot(121); hold on; axis square;
        var.handles(subject).extraLikeTitle = title('Likelihood h(c)');
        var.handles(subject).extraLikeXLab = xlabel('Contrast');
        var.handles(subject).extraPriorAx = subplot(122); hold on; axis square;
        var.handles(subject).extraPriorTitle = title('Prior');
        var.handles(subject).extraPriorXLab = xlabel('Speed (deg s^{-1})');
        
        % Pretty up the axes
        set(var.handles(subject).extraLikeAx, 'XLim', [.04 1], 'YLim', [.5 3], 'XScale', 'log', ...
            'TickDir', 'out', 'FontSize', 12, 'FontName', 'Helvetica');
        set(var.handles(subject).extraPriorAx, 'XLim', [.4 13], 'YLim', [1e-14 2], 'XScale', 'log', ...
            'YScale', 'log', 'TickDir', 'out', 'FontSize', 12, 'FontName', 'Helvetica');
        
        % Plot each resampled likelihood
        % First likelihood on real data
        set(var.handles(subject).extraFig, 'CurrentAxes', var.handles(subject).extraLikeAx);
        var.handles(subject).extraLikeLine = ...
        plot(var.data(subject).testContrasts, var.params(subject).likeWidth, ...
        'LineStyle', '--', 'LineWidth', 2, 'Marker', 'none', 'Color', [1 .5 .5]);
    
        % Then resampled data
        if var.flags.fitG
        else
            var.handles(subject).extraLikeLines = ...
                plot(repmat(var.data(subject).testContrasts, 1, size(var.params(subject).hHat, 1)), ...
                var.params(subject).hHat' ./ var.data(subject).G, ...
                'LineStyle', '--', 'Marker', 'none', 'Color', [.75 .75 .75]);
        end
        
        % Plot each resampled prior
        % First prior on real data
        set(var.handles(subject).extraFig, 'CurrentAxes', var.handles(subject).extraPriorAx);
        var.handles(subject).extraPriorLine = ...
        plot(var.data(subject).refVels, var.params(subject).prior(1, :), ...
        'LineStyle', '--', 'Marker', 'none', 'Color', [1 .5 .5]);
    
        % Then resampled data
        var.handles(subject).extraPriorLines = ...
            plot(repmat(var.data(subject).refVels, 1, size(var.params(subject).prior, 1) - 1), ...
            var.params(subject).prior(2:end, :)', ...
            'LineStyle', '--', 'Marker', 'none', 'Color', [.75 .75 .75]);

    end
end

%% Organize
var = organizeVar(var);