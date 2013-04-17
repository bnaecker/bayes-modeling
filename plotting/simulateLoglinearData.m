%% simulate the observer
% define stuff
nTrials = 50000;
testVels = [.5 1 2 4 8 12];
refVels = testVels;
testContrasts = [.05 .1 .2 .4 .5];
refContrasts = [.05 .5];
maxSigma = .1;

%% make reference velocities and contrasts
c1inds = randi(length(refContrasts), nTrials, 1);
c2inds = randi(length(testContrasts), nTrials, 1);
c1 = refContrasts(c1inds)';
c2 = testContrasts(c2inds)';
sig1 = (1 ./ c1) .* maxSigma;
sig2 = (1 ./ c2) .* maxSigma;
v1inds = randi(length(refVels), nTrials, 1);
v1 = refVels(v1inds)';

%%% OLD %%%%
% v2inds = randi(length(testVels), nTrials, 1);
% v2 = testVels(v2inds)';

%% Make test elocities
% test velocities should be distributed around the reference velocities to
% approximate the distribution given by a staircase procedure. But the
% distribution should not be too dense. So pick a handful of points around
% the velocities, and sample with a broadly gaussian distribution.
minFrac = .5; maxFrac = 1.5;
nFracs = 5;
f = linspace(minFrac, maxFrac, nFracs)';
f = f * ones(1, length(testVels));
testV = ones(nFracs, 1) * testVels;
testV = testV .* f;

% Sample from testV appropriately
sigma = 1;
weight = normpdf(linspace(-.5, .5, nFracs), 0, sigma);
v2 = zeros(size(v1));
for ti = 1:nTrials
    vs = testV(:, v1inds(ti));
    v2(ti) = randsample(vs, 1, true, weight);
end

%% Add gaussian noise to the samples
n1 = randn(nTrials, 1) .* sig1;
n2 = randn(nTrials, 1) .* sig2;
vel1 = v1 + n1;
vel2 = v2 + n2;

%% prior description
priorSlopes = logspace(-.5, -1, length(testVels));
prior1 = priorSlopes(v1inds)';
% prior2 = priorSlopes(v2inds)';
prior2 = interp1(testVels, priorSlopes, vel2, 'linear', 'extrap');

%% compute posterior estimates
alpha1 = prior1 .* sig1 .^ 2;
alpha2 = prior2 .* sig2 .^ 2;
p1 = vel1 - alpha1;
p2 = vel2 - alpha2;

%% choices are maxima of these
T = p2 > p1;

%% save choices to run through fitting code
s4 = [v1 c1 v2 c2 nan(nTrials, 4) T nan(nTrials, 2)];
save('/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/NewCode/Data/s4.mat', ...
    's4', 'sig1', 'sig2');

% save information about distributions, for plotting/comparisons
contrasts = testContrasts;
sigmas = (1 ./ testContrasts) .* maxSigma;
vels = testVels;
save(['/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/NewCode/Data/' ...
    'loglinearDists.mat'], 'sigmas', 'vels', 'priorSlopes', 'contrasts', 'sig1', 'sig2');