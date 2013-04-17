function ds = plotFittedDistributions(ds)
%
% FUNCTION ds = plotFittedDistributions(ds)
%
% Plots the actual distributions fitted to the Bayesian ideal observer
% model, both on a linear and logarithmic scale.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

%% notify
fprintf('Plotting fitted distributions ... ');

%% loop over subjects
for si = 1:ds.info.nSubjects
    % preallocate handles to lines
    ds.handles(si).priorDists = zeros(ds.info.nUniqueContrasts, 2);
    ds.handles(si).likeDists = zeros(ds.info.nUniqueRefVels, ...
        ds.info.nUniqueContrasts, 2);
    ds.handles(si).postDists = ds.handles(si).likeDists;
    ds.handles(si).distAx = zeros(ds.info.nUniqueContrasts, 2);
    
    % make a figure
    ds.handles(si).linearDistFig = figure;
    if ds.info.isSimData
        ext = ', simulated data';
    else
        ext = '';
    end
    set(ds.handles(si).linearDistFig, 'Color', 'w', 'NumberTitle', 'off', ...
        'Name', ['prior and likelihood distributions, linear scale ' ...
        'subject ' num2str(si) ext]);
    
    % transform speed back to linear axis
    if ds.info.isSimData
        ax = ds.params(si).interpAx;
        vels = ds.data(si).refVels;
    else
        ax = ds.velTrans.iTransFun(ds.params(si).interpAx);
        vels = ds.velTrans.iTransFun(ds.data(si).refVels);
    end
    mx = 15; % max deg/s to plot
    plotInds = ax < mx;
    ax = ax(plotInds);
    
    % get spacing
    dx = diff(ax(1:2));
    
    % determine good number of subplots
    nSubplots = numSubplots(ds.info.nUniqueContrasts);
    
    % plot fits on the linear axis
    for ci = 1:ds.info.nUniqueContrasts
        ds.handles(si).distAx(ci, 1) = subplot(nSubplots(1), nSubplots(2), ci);
        hold on; axis square;
        
        % plot prior
        ds.handles(si).priorDists(ci, 1) = ...
            plot(ax, ds.params(si).interpPrior(plotInds, 1), ...
            'LineWidth', 2, 'Color', 'k', 'Marker', 'none');
        
        % plot each likelihood and posterior
        for vi = 1:ds.info.nUniqueRefVels
            % likelihood
            ll = normpdf(ax, vels(vi), ds.params(si).likeWidth(ci, 1))';
            ll = ll ./ (sum(ll) * dx);
            ds.handles(si).likeDists(vi, ci, 1) = plot(ax, ll, ...
                'LineWidth', 1.5, 'Marker', 'none', 'Color', [0.8 0.1 0.1]);
            
            % posterior
            post = ds.params(si).interpPrior(plotInds, 1) .* ll;
            post = post ./ (sum(post) * dx);
            ds.handles(si).postDists(vi, ci, 1) = plot(ax, post(plotInds), ...
                'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.2 0.2 0.6]);
        end
        
        % pretty up the axes
        set(ds.handles(si).distAx(ci, 1), 'TickDir', 'out', 'TickLength', [.003 .003], ...
            'FontSize', 12, 'FontWeight', 'bold', 'XLim', [1e-1 mx +  0.2 * mx]);
        
        % labels and title
        ds.handles(si).distAxXLabel(ci, 1) = xlabel('speed (deg s^{-1})', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).distAxYLabel(ci, 1) = ylabel('probability', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).distAxTitle(ci, 1) = title(...
            ['contrast = ' num2str(ds.data(si).testContrasts(ci), '%.2f')], ...
            'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % make a figure
    ds.handles(si).logDistFig = figure;
    set(ds.handles(si).logDistFig, 'Color', 'w', 'NumberTitle', 'off', ...
        'Name', ['prior and likelihood distributions, logarithmic scale ' ...
        'subject ' num2str(si) ext]);
    
    % plot the fits on the log axis
    for ci = 1:ds.info.nUniqueContrasts
        ds.handles(si).distAx(ci, 2) = subplot(nSubplots(1), nSubplots(2), ci);
        hold on; axis square;
        
        % plot prior
        ds.handles(si).priorDists(ci, 2) = ...
            plot(ax, ds.params(si).interpPrior(plotInds, 1), ...
            'LineWidth', 2, 'Color', 'k', 'Marker', 'none');
        
        % plot each likelihood and posterior
        for vi = 1:ds.info.nUniqueRefVels
            % likelihood
            ll = normpdf(ax, vels(vi), ds.params(si).likeWidth(ci, 1))';
            ll = ll ./ (sum(ll) * dx);
            ds.handles(si).likeDists(vi, ci, 2) = plot(ax, ll, ...
                'LineWidth', 1, 'Marker', 'none', 'Color', [0.8 0.1 0.1]);
            
            % posterior
            post = ds.params(si).interpPrior(plotInds, 1) .* ll;
            post = post ./ (sum(post) * dx);
            ds.handles(si).postDists(vi, ci, 2) = plot(ax, post(plotInds), ...
                'LineWidth', 1, 'LineStyle', '--', 'Color', [0.2 0.2 0.6]);
        end
        
        % pretty up the axes
        set(ds.handles(si).distAx(ci, 2), 'TickDir', 'out', 'TickLength', [.003 .003], ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'XScale', 'log', 'YScale', 'log', ...
            'XLim', [1e-1 mx + 0.2 * mx], 'YLim', [1e-20 1e2]);
        
        % labels and title
        ds.handles(si).distAxXLabel(ci, 2) = xlabel('speed (deg s^{-1})', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).distAxYLabel(ci, 2) = ylabel('probability', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).distAxTitle(ci, 2) = title(...
            ['contrast = ' num2str(ds.data(si).testContrasts(ci), '%.2f')], ...
            'FontSize', 14, 'FontWeight', 'bold');
    end
end

%% order
ds = orderDataStruct(ds);
fprintf('done.\n');