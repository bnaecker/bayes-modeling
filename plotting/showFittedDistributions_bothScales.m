function data = showFittedDistributions_bothScales(data)
%
% data = showFittedDistributions_bothScales(data)
%
% Helper function to show the fitted prior and likelihoods on the same axes
% with both log- and linearly-scaled y-axes.
%
% See also: RunSpeedDiscriminationModel.m, showSimDists.m,
% showFittedDistributions.m
%
% (c) bnaecker@stanford.edu 4 Apr 2012

%% Parse input
if ~isstruct(data)
    error('showFittedDistributions_bothScales:badInput', ...
        'Input must be a data structure returned by RunSpeedDiscriminationModel.m');
end

%% Setup plot colors
priorColor = 'k';
likeColor = [.8 .1 .1];

%% Plot
nSubplots = [data.info.nUniqueContrasts 2];
for si = 1:data.info.nSubjects
    % setup handles
    data.handles(si).allPriors = zeros(data.info.nBoots, 2);
    data.handles(si).allLikelihoods = zeros(data.info.nUniqueRefVels, ...
        data.info.nUniqueContrasts, 2);
    data.handles(si).allPosteriors = data.handles(si).allLikelihoods;
    data.handles(si).showAx = zeros(data.info.nUniqueContrasts, 2);
    
    % make figure
    data.handles(si).distFig = figure;
    set(data.handles(si).distFig, 'Color', 'w', 'NumberTitle', 'off', ...
        'Name', ['Priors and likelihoods, Subject ' num2str(si)]);
    
    % transform the speed axis back to linear
    if data.info.isSimData
        ax = data.params(si).interpAx;
        vels = data.data(si).refVels;
    else
        ax = data.vTrans.iTransFun(data.params(si).interpAx);
        vels = data.vTrans.iTransFun(data.data(si).refVels);
    end
    dd = diff(ax); dx = dd(1);
    
    % plot the fits on the linear axis (left column)
    for ci = 1:data.info.nUniqueContrasts
        data.handles(si).showAx(ci, 1) = subplot(nSubplots(1), 2, 2 * ci - 1);
        hold on;
        % Prior
        data.handles(si).allPriors(1, 1) = ...
                plot(ax, data.params(si).interpPrior(1, :), ...
                'LineWidth', 2, 'Marker', 'none', 'Color', priorColor);
        for vi = 1:data.info.nUniqueRefVels
            % likelihood
            ll = normpdf(ax, vels(vi), data.params(si).likeWidth(1, ci));
            ll = ll ./ (sum(ll) * dx);
            data.handles(si).allLikelihoods(vi, ci, 1) = plot(ax, ll, ...
                'LineWidth', 2, 'Marker', 'none', 'Color', likeColor);
            
            % posterior
            post = data.params(si).interpPrior(1, :) .* ll;
            post = post ./ (sum(post) * dx);
            data.handles(si).allPosteriors(vi, ci, 1) = plot(ax, post, ...
                'LineWidth', 1', 'Marker', 'none', 'Color', 'g');
        end
        
        % pretty up axes
        set(data.handles(si).showAx(ci, 1), 'TickDir', 'out', 'TickLength', [.005 0], ...
            'XScale', 'linear', 'YScale', 'linear');
%         data.handles(si).showAxTitle(1) = title(['Contrast = ' ...
%             num2str(data.data(1).testContrasts(ci))]);
        if ci == 1
            title('Linear');
        end
        data.handles(si).showAxYLabel(1) = ylabel(['Contrast = ' ...
            num2str(data.data(1).testContrasts(ci))]);
        if ci == data.info.nUniqueContrasts
            data.handles(si).showAxXLabel(ci) = xlabel('Speed (deg/s)');
        end
    end
    
    % plot the fits on the log axis
    for ci = 1:data.info.nUniqueContrasts
        data.handles(si).showAx(ci, 2) = subplot(nSubplots(1), 2, 2 * ci);
        hold on;
        % Prior
        data.handles(si).allPriors(1, 2) = ...
                plot(ax, data.params(si).interpPrior(1, :), ...
                'LineWidth', 2, 'Marker', 'none', 'Color', priorColor);
        for vi = 1:data.info.nUniqueRefVels
            % likelihood
            ll = normpdf(ax, vels(vi), data.params(si).likeWidth(1, ci));
            ll = ll ./ (sum(ll) * dx);
            data.handles(si).allLikelihoods(vi, ci, 2) = plot(ax, ll, ...
                'LineWidth', 2, 'Marker', 'none', 'Color', likeColor);
            
            % posterior
            post = data.params(si).interpPrior(1, :) .* ll;
            post = post ./ (sum(post) * dx);
            data.handles(si).allPosteriors(vi, ci, 2) = plot(ax, post, ...
                'LineWidth', 1', 'Marker', 'none', 'Color', 'g');
        end
        
        % pretty up axes
        set(data.handles(si).showAx(ci, 2), 'TickDir', 'out', 'TickLength', [.005 0], ...
            'XScale', 'linear', 'YScale', 'log', 'YLim', [1e-10 1e2]);
%         data.handles(si).showAxTitle(2) = title(['Contrast = ' ...
%             num2str(data.data(1).testContrasts(ci))]);
        if ci == 1
            title('Log');
        end
        if ci == data.info.nUniqueContrasts
            data.handles(si).showAxXLabel(ci) = xlabel('Speed (deg/s)');
        end
    end
    
    % quick check that we're in a reasonable range
    lims = get(data.handles(si).showAx, 'XLim');
    limit = 20;
    for ii = 1:length(lims)
        lims{ii}( lims{ii} > limit ) = limit;
    end
    set(data.handles(si).showAx, {'XLim'}, lims);
end