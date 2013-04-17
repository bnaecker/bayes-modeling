function data = showFittedDistributions(data, varargin)
%
% function data = showFittedDistributions(data)
%
% Helper function to plot the fitted prior and likelihoods on the same
% axis, to visualize all fits and to compare the two distributions
%
% See also: RunSpeedDiscriminationModel.m, showSimDists.m
%
% (c) bnaecker@stanford.edu 7 Mar 2012

%% Parse input
if nargin == 2 && islogical(varargin{1})
    plotResamples = varargin{1};
    axScale = 'linear';
elseif nargin == 2 && ischar(varargin{1})
    plotResamples = false;
    axScale = varargin{1};
else
    plotResamples = false;
    axScale = 'linear';
end

%% Setup plot colors
priorColor = 'k';
likeColor = [.8 .1 .1];

%% Plot
nSubplots = numSubplots(data.info.nUniqueContrasts);
for si = 1:data.info.nSubjects
    % setup handles
    data.handles(si).allPriors = zeros(data.info.nBoots, 1);
    data.handles(si).allLikelihoods = zeros(data.info.nUniqueRefVels, ...
        data.info.nUniqueContrasts);
    data.handles(si).allPosteriors = data.handles(si).allLikelihoods;
    data.handles(si).showAx = zeros(data.info.nUniqueContrasts, 1);
    
    % make a figure
    data.handles(si).distFig = figure; hold on;
    set(data.handles(si).distFig, 'Color', 'w', 'NumberTitle', 'off', ...
        'Name', ['Priors and likelihood widths, Subject ' num2str(si)]);
    
    % transform the speed axis back to linear
    if data.info.isSimData
        ax = data.params(si).interpAx;
        vels = data.data(si).refVels;
    else
        ax = data.vTrans.iTransFun(data.params(si).interpAx);
        vels = data.vTrans.iTransFun(data.data(si).refVels);
    end
    dd = diff(ax); dx = dd(1);

    % plot the fits to original data with thicker lines
    for ci = 1:data.info.nUniqueContrasts
        data.handles(si).showAx(ci) = subplot(nSubplots(2), nSubplots(1), ci);
        hold on;
        % Prior
        data.handles(si).allPriors(1) = ...
                plot(ax, data.params(si).interpPrior(1, :), ...
                'LineWidth', 2, 'Marker', 'none', 'Color', priorColor);
        for vi = 1:data.info.nUniqueRefVels
            % likelihood
            ll = normpdf(ax, vels(vi), data.params(si).likeWidth(1, ci));
            ll = ll ./ (sum(ll) * dx);
            data.handles(si).allLikelihoods(vi, ci) = plot(ax, ll, ...
                'LineWidth', 2, 'Marker', 'none', 'Color', likeColor);
            
            % posterior
            post = data.params(si).interpPrior(1, :) .* ll;
            post = post ./ (sum(post) * dx);
            data.handles(si).allPosteriors(vi, ci) = plot(ax, post, ...
                'LineWidth', 1', 'Marker', 'none', 'Color', 'g');
        end
        
        % pretty up axes
        set(data.handles(si).showAx(ci), 'TickDir', 'out', 'TickLength', [.005 0], ...
            'XScale', 'linear', 'YScale', axScale);
        data.handles(si).showAxTitle = title(['Contrast = ' ...
            num2str(data.data(1).testContrasts(ci))]);
        data.handles(si).showAxXLabel = xlabel('Speed (deg/s)');
        data.handles(si).showAxYLabel = ylabel('pdf');
    end
    
    % plot the fits to resampled data
    if plotResamples
        data.handles(si).allPriors(2:end) = ...
            plot(ax, data.params(si).interpPrior(2:end, :)', ...
            'LineWidth', 1, 'Marker', 'none', 'Color', priorColor);
        for vi = 1:data.info.nUniqueRefVels
            for bi = 2:data.info.nBoots
                data.handles(si).allLikelihoods(bi, vi) = ...
                    plot(ax, normpdf(ax, vels(vi), ...
                    data.params(si).likeWidth(bi, vi)), ...
                    'LineWidth', 1, 'Marker', 'none', 'Color', likeColor);
            end
        end
    end
end