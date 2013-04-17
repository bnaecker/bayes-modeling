function [sim ds] = testModelOnSimulatedData(priorType)
%
% FUNCTION [sim ds] = testModelOnSimulatedData(priorType)
%
% This function generates fake 2AFC choice data from a Bayesian ideal
% observer assuming the requested prior type, then fits the model using the
% including fitting code, and compares the results.
%
% (c) bnaecker@stanford.edu 13 Nov 2012

%% simulate the data
sim = simulateDiscriminationData('priorType', priorType);

%% fit the model to the simulated data
ds = runSpeedDiscriminationModel('nSubjects', 0, 'priorType', priorType);

%% generate plots overlaying the fitted distributions
% prior
figure, set(gcf, 'Color', 'w', 'NumberTitle', 'off', ...
    'Name', 'Comparison of prior distribution');
if strcmp(priorType, 'gaussian')
    plot(ds.params.interpAx, normpdf(ds.params.interpAx, 0, ds.params.gamma) ./ ...
        (sum(normpdf(ds.params.interpAx, 0, ds.params.gamma)) * 0.01), 'Color', 'k', ...
        'LineWidth', 2); hold on;
    plot(ds.params.interpAx, normpdf(ds.params.interpAx, 0, sim.gamma) ./ ...
        (sum(normpdf(ds.params.interpAx, 0, sim.gamma)) * 0.01), ...
        'Color', 'r', 'LineWidth', 2);
elseif strcmp(priorType, 'loglinear');
    plot(ds.params.interpAx, ds.params.interpPrior, 'Color', 'k', ...
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
plot(ds.data.testContrasts, ds.params.likeWidth, ...
    'Color', 'k', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 6); hold on;
plot(sim.contrasts, sim.likeWidth, ...
    'Color', 'r', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 6);
legend({'fitted', 'true'}, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('contrast', 'FontSize', 12, 'FontWeight', 'bold');
box off; grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.003 0.003], ...
    'FontSize', 12, 'FontWeight', 'bold');