function h = showSimDists(varargin)
%
% function h = showSimDists(varargin)
%
% Helper function to plot the distributions from which simulated ideal
% observer samples were drawn (likelihoods).
%
% See also: simulateGaussianData.m, showFittedDistributions.m,
% RunSpeedDiscriminationModel.m
%
% (c) bnaecker@stanford.edu 26 Mar 2012

%% parse input
baseDir = ['/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination' ...
    '/NewCode/Data'];
if nargin == 0
    warning('showSimDists:noInput', ...
        'No simulation type given, assuming "gaussian"');
    priorType = 'gaussian';
elseif nargin == 1 && ...
        (strcmp(varargin{1}, 'gaussian') || strcmp(varargin{1}, 'loglinear'))
    priorType = varargin{1};
elseif nargin == 2 && ...
        (strcmp(varargin{1}, 'gaussian') || strcmp(varargin{1}, 'loglinear'))
    priorType = varargin{1};
    axScale = varargin{2};
else
    % error for now
    error('showSimDists:uknownPriorType', ...
        ['The requested prior type is unknown. As of %s, can only accept' ...
        '"loglinear" or "gaussian" as arguments'], datestr(now, 'dd mmm yyyy'));
end

%% Gather data on distributions, axes
fileName = fullfile(baseDir, [priorType 'Dists.mat']);
s = load(fileName);
sigmas = s.sigmas;
vels = s.vels;
contrasts = s.contrasts;
dx = .01;
speedAx = dx:dx:(max(vels) + .5 * max(vels));
if strcmp(priorType, 'gaussian')
    priorDef = s.gamma;
else
    priorDef = s.priorSlopes;
end

%% setup handle output structure
h = struct('handles', ...
           struct('figure', 0, ...
                  'axes', zeros(length(sigmas), 1), ...
                  'priorLine', zeros(length(sigmas), 1), ...
                  'likeLines', zeros(length(sigmas), length(vels)), ...
                  'xlabels', zeros(length(sigmas), 1), ...
                  'ylabels', zeros(length(sigmas), 1), ...
                  'titles', zeros(length(sigmas), 1)), ...
           'simData', ...
           struct('sigmas', sigmas, ...
                  'vels', vels, ...
                  'priorDef', priorDef, ...
                  'contrasts', contrasts', ...
                  'prior', zeros(size(speedAx)), ...
                  'likelihoods', {cell(length(sigmas), length(vels))}));
       
%% Plot
nPlots = numSubplots(length(contrasts));
h.handles.figure = figure;
set(h.handles.figure, 'Color', 'w', 'Name', ...
    ['Simulated prior and likelihoods, ' priorType], 'NumberTitle', 'off');

% one subplot for each contrast
for si = 1:length(contrasts)
    % make axes, labels
    h.handles.axes(si) = subplot(nPlots(2), nPlots(1), si); hold on;
    h.handles.xlabels(si) = xlabel('Speed (deg s^{-1})');
    h.handles.ylabels(si) = ylabel('pdf');
    h.handles.titles(si) = title(['Contrast = ' num2str(contrasts(si))]);
    
    % plot prior
    if strcmp(priorType, 'gaussian')
        pr = normpdf(speedAx, 0, priorDef);
        pr = pr ./ (sum(pr) * dx);
    else
        slopes = interp1(vels, priorDef, speedAx, 'linear', 'extrap');
        prCum = cumsum(slopes) * dx;
        pr = exp(-prCum);
        pr = pr ./ (sum(pr) * dx);
    end
    h.handles.priorLine(si) = plot(speedAx, pr, ...
        'Color', 'k', 'LineWidth', 2, 'Marker', 'none');
    % save on first go-round
    if si == 1
        h.simData.prior = pr;
    end
    
    % plot likelihood for this contrast at each velocity
    for vi = 1:length(vels)
        ll = normpdf(speedAx, vels(vi), sigmas(si));
        ll = ll ./ (sum(ll) * dx);
        h.simData.likelihoods{si, vi} = ll;
        h.handles.likeLines(si, vi) = plot(speedAx, ll, ...
            'Color', [.8 .1 .1], 'Marker', 'none', 'LineWidth', 2);
    end
    
    % pretty up axes
    if nargin < 2
        axScale = 'linear';
    end
    set(h.handles.axes(si), 'TickDir', 'out', 'TickLength', [.005 0], ...
        'XScale', 'linear', 'YScale', 'linear'); box off;
end
h.simData.speedAx = speedAx;