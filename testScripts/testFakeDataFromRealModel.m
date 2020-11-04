function [sim ds] = testFakeDataFromRealModel(priorType)
%
% FUNCTION [sim ds] = testFakeDataFromRealModel(varargin)
%
% This function fits the Bayesian ideal observer model to real data, then
% generates fake data with the fitted characteristics, and then fits the
% model again to that simulated data.
%
% (c) bnaecker@stanford.edu 14 Nov 2012

%% fit model
fprintf('\nFitting model to real data ... ');
ds = runSpeedDiscriminationModel('nSubjects', 1, 'priorType', priorType);

%% simulate
if strcmp(priorType, 'gaussian')
    priorDefStr = 'gamma';
    priorDef = ds.params.gamma(1);
else
    priorDefStr = 'slopes';
    priorDef = ds.params.slopeHat(:, 1);
end
sim = simulateDiscriminationData('priorType', priorType, ...
    priorDefStr, priorDef);

%% test on data simulated with parameters from the fitted model
fprintf('\nFitting model to data simulated from real fitted model ... ');
ds = runSpeedDiscriminationModel('nSubjects', 0, 'priorType', priorType);

%% generate plots overlaying the fitted distributions
% prior
figure, set(gcf, 'Color', 'w', 'NumberTitle', 'off', ...
    'Name', 'Comparison of prior distribution');
if strcmp(sim.priorType, 'gaussian')
    plot(ds.params.interpAx, normpdf(ds.params.interpAx, zeros(size(ds.params.gamma)), ...
		ds.params.gamma) ./ ...
        (sum(normpdf(ds.params.interpAx, zeros(size(ds.params.gamm)), ...
		ds.params.gamma)) * 0.01), 'Color', 'k', ...
        'LineWidth', 2); hold on;
    plot(ds.params.interpAx, normpdf(ds.params.interpAx, zeros(size(sim.gamma)), sim.gamma) ./ ...
        (sum(normpdf(ds.params.interpAx, zeros(size(sim.gamma)), sim.gamma)) * 0.01), ...
        'Color', 'r', 'LineWidth', 2);
elseif strcmp(sim.priorType, 'loglinear');
    plot(ds.params.interpAx, ds.params.interpPrior(:, 1), 'Color', 'k', ...
        'LineWidth', 2); hold on;
    plot(sim.vAx, sim.prior, 'Color', 'r', 'LineWidth', 2);
end
legend({'fitted', 'true'}, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('speed (deg s^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
box off; grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.003 0.003], ...
    'FontSize', 12, 'FontWeight', 'bold');

% likelihoods
figure, set(gcf, 'Color', 'w', 'NumberTitle', 'off', ...
    'Name', 'Comparison of likelihood distributions');
plot(ds.data.testContrasts, ds.params.likeWidth(:, 1), ...
    'Color', 'k', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 6); hold on;
plot(sim.contrasts, sim.likeWidth, ...
    'Color', 'r', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 6);
legend({'fitted', 'true'}, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('contrast', 'FontSize', 12, 'FontWeight', 'bold');
box off; grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.003 0.003], ...
    'FontSize', 12, 'FontWeight', 'bold');
