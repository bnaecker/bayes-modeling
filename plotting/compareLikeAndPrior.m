function hands = compareLikeAndPrior(data, h)
%
% hands = compareLikeAndPrior(data, h)
%
% Helper function to show the simulated and fitted likelihood widths and
% prior slopes/widths
%
% (c) bnaecker@stanford.edu 3 Apr 2012


%% Setup fig
hands.figure = figure;
set(hands.figure, 'Color', 'w', 'NumberTitle', 'off', ...
    'Name', 'Comparison of simulated and fitted distributions');

%% compute error
likeErr = sum((h.simData.sigmas - data.params.likeWidth(1, :)) .^ 2);
dx = .01;
simAx = max(min(h.simData.speedAx), min(data.params.interpAx)) : dx : ...
    min(max(h.simData.speedAx), max(data.params.interpAx));
interpSimPrior = interp1(h.simData.speedAx, h.simData.prior, simAx, 'linear');
interpFitPrior = interp1(data.params.interpAx, data.params.interpPrior(1, :), ...
    simAx, 'linear');
prErr = sum((interpSimPrior - interpFitPrior) .^ 2);

%% plot likelihood widths
hands.ax(1) = subplot(121);
xx = data.data.testContrasts;
hands.lines(1:2) = plot(xx, data.params.likeWidth(1, :), 'r', ...
    xx, h.simData.sigmas, 'k');
hands.legend(1) = legend('fitted', 'true');
hands.title(1) = title('Likelihood widths');
hands.xlab(1) = xlabel('Contrast'); 
set(hands.ax(1), 'TickDir', 'out'); box off;
hands.text(1) = text(.5, 1, ['SSE = ' num2str(likeErr)]);

%% plot prior
hands.ax(2) = subplot(122);
hands.lines(3:4) = plot(data.params.interpAx, data.params.interpPrior(1, :), 'r', ...
    h.simData.speedAx, h.simData.prior, 'k');
hands.legend(1) = legend('fitted', 'true');
hands.title(1) = title('Priors');
hands.xlab(1) = xlabel('Speed (deg s^{-1})'); 
set(hands.ax(2), 'TickDir', 'out'); box off;
hands.text(2) = text(10, .14, ['SSE \approx' num2str(prErr)]);