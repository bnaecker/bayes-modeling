function ds = plotPriorFit(ds)
%
% FUNCTION ds = plotPriorFit(ds)
%
% Plots the actual prior distribution and the widths of the likelihood
% function for each subject. Corresponds to figure 4 in Stocker &
% Simoncelli 2006.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

% notify
fprintf('Plotting prior and likelihood widths ... ');

% loop over subjects
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    % make figure and supblot for likelihood widths
    ds.handles(si).priorFig = figure;
    if ds.info.isSimData
        ext = ', simulated data';
    else
        ext = '';
    end
    set(ds.handles(si).priorFig, 'Color', 'w', 'NumberTitle', 'off', ...
        'Name', ['Subject ' num2str(si) ext]);
    ds.handles(si).likeWidthAx = subplot(121); hold on; axis square;
    ds.handles(si).likeWidthTitle = title('width of likelihood', ...
        'FontSize', 14, 'FontWeight', 'bold');
    ds.handles(si).likeWidthXLabel = xlabel('contrast', ...
        'FontSize', 12, 'FontWeight', 'bold');
    set(ds.handles(si).likeWidthAx, 'TickDir', 'out', ...
        'TickLength', [.003 .003], 'FontSize', 12, 'FontWeight', 'bold', ...
        'XScale', 'log');
    
    
    % plot likelihood and error bars, if requested
    ds.handles(si).likeWidthLine = ...
        plot(ds.data(si).testContrasts, ds.params(si).likeWidth(:, 1), ...
        'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'k');
    if ds.flags.plotSmears
        ds.handles(si).likeWidthSmears = ...
            fill([ds.data(si).testContrasts; flipud(ds.data(si).testContrasts)], ...
            [ds.stats(si).likeWidthLB; flipud(ds.stats(si).likeWidthUB)], ...
            0.75 .* [1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    % plot functional form to likelihood, if computed
    if ds.flags.fitLikeFun
        if ds.flags.useCFit
            ds.handles(si).likeFunLine = ...
                plot(ds.likeOpt(si).domain, ...
                ds.likeOpt(si).fitObj(ds.likeOpt(si).domain), ...
                'LineStyle', '--', 'Marker', 'none', 'Color', 'k');
        else
            ds.handles(si).likeFunLine = ...
                plot(ds.likeOpt(si).domain, ...
                ds.likeOpt(si).likeFunVals, ...
                'LineStyle', '--', 'Marker', 'none', 'Color', 'k');
        end
    end
    
    
    % subplot for prior
    ds.handles(si).priorAx = subplot(122); hold on; axis square;
    ds.handles(si).priorTitle = title('prior', ...
        'FontSize', 12, 'FontWeight', 'bold');
    ds.handles(si).priorXLabel = xlabel('speed (deg s^{-1})', ...
        'FontSize', 12, 'FontWeight', 'bold');
    set(ds.handles(si).priorAx, 'TickDir', 'out', ...
        'TickLength', [.003 .003], 'FontSize', 12, 'FontWeight', 'bold', ...
        'XScale', 'log', 'YScale', 'log');
    
    % plot prior and error bars, if requested
    if ds.info.isSimData
        vAx = ds.data(si).refVels;
    else
        vAx = ds.velTrans.iTransFun(ds.data(si).refVels);
    end
    ds.handles(si).priorLine = ...
        plot(vAx, ...
        ds.params(si).prior(:, 1), ... 
        'Color', 'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6);
    if ds.flags.plotSmears
        ds.handles(si).priorSmears = ...
            fill([vAx; flipud(vAx)], ...
            [ds.stats(si).priorLB(:, 1); flipud(ds.stats(si).priorUB(:, 1))], ...
            0.75 .* [1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    
    % plot individual resamples if requested
    if ds.flags.plotStraps
        % create figure, axes, and titles/labels
        ds.handles(si).resampleFig = figure;
        set(ds.handles(si).resampleFig, 'Color', 'w', 'NumberTitle', 'off', ...
            'Name', ['Subject ' num2str(si), ', fits to resampled data']);
        ds.handles(si).resampleLikeAx = subplot(121); hold on;
        ds.handles(si).resampleLikeTitle = title('width of likelihood', ...
            'FontSize', 14, 'FontWeight', 'bold');
        ds.handles(si).resampleLikeXlabel = xlabel('contrast', ...
            'FontSize', 12, 'FontWeight', 'bold');
        set(ds.handles(si).resampleLikeAx, 'TickDir', 'out', 'TickLength', [.003 .003], ...
            'FontSize', 12, 'FontWeight', 'bold');
        
        % plot likelihood resamples
        % first to original data, in red
        ds.handles(si).resampledLikeLine = ...
            plot(ds.data(si).testContrasts, ds.params(si).likeWidth(:, 1), ...
            'Color', [1 .5 .5], 'LineStyle', '--', 'Marker', 'none', 'LineWidth', 2);
        % individual resamples, in black
        if ds.flags.fitLikeSpeed
            % optional functionality
        else
            ds.handles(si).resampledLikeLines = ...
                plot(ds.data(si).testContrasts * ones(1, size(ds.params(si).likeWidth, 2) - 1), ...
                ds.params(si).likeWidth(:, 2:end), ...
                'LineStyle', '--', 'LineWidth', 1, 'Marker', 'none', 'Color', 0.75 .* [1 1 1]);
        end
        
        % then prior
        ds.handles(si).resamplePriorAx = subplot(122); hold on;
        ds.handles(si).resamplePriorTitle = title('prior', ...
            'FontSize', 14, 'FontWeight', 'bold');
        ds.handles(si).resamplePriorXlabel = xlabel('speed (deg s^{-1})', ...
            'FontSize', 12, 'FontWeight', 'bold');
        set(ds.handles(si).resamplePriorAx, 'TickDir', 'out', 'TickLength', [.003 .003], ...
            'FontSize', 12, 'FontWeight', 'bold', 'YScale', 'log', 'XScale', 'log');
        
        % plot prior resamples
        % first to original data, in red
        ds.handles(si).resampledPriorLine = ...
            plot(vAx, ds.params(si).prior(:, 1), ...
            'Color', [1 .5 .5], 'LineStyle', '--', 'Marker', 'none', 'LineWidth', 2);
        % individual resamples, in black
        ds.handles(si).resampledPriorLines = ...
            plot(vAx * ones(1, size(ds.params(si).prior, 2) - 1), ...
            ds.params(si).prior(:, 2:end), ...
            'LineStyle', '--', 'LineWidth', 1, 'Marker', 'none', 'Color', 0.75 .* [1 1 1]);
    end
end

%% order
ds = orderDataStruct(ds);
fprintf('done.\n');