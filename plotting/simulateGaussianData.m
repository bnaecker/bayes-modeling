%% simulate the observer
% define stuff
nTrials = 50000;
testVels = [.5 1 2 4 8 12];
refVels = testVels;
testContrasts = [.05 .1 .2 .4 .5];
refContrasts = [.05 .5];
% testContrasts = refContrasts;
gamma = 5; % prior sd
maxSigma = .5;   % max likelihood width

%% make reference velocities and contrasts
c1inds = randi(length(refContrasts), nTrials, 1);
c2inds = randi(length(testContrasts), nTrials, 1);
c1 = refContrasts(c1inds)';
c2 = testContrasts(c2inds)';
sig1 = (1 ./ c1) .* maxSigma;
sig2 = (1 ./ c2) .* maxSigma;
% sigs = [3 6];
% sig1 = sigs(c1inds)';
% sig2 = sigs(c2inds)';
v1 = refVels(randi(length(refVels), nTrials, 1))';

%% Make test elocities

%%%%% REVISIT THIS. IT MAY NOT BE THE BEST WAY TO GENERATE DATA FOR THE
%%%%% GAUSSIAN CASE



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

% vel1 = max(v1 + n1, 0);
% vel2 = max(v2 + n2, 0);

%% Add gaussian noise to the samples
n1 = randn(nTrials, 1) .* sig1;
n2 = randn(nTrials, 1) .* sig2;
vel1 = v1 + n1;
vel2 = v2 + n2;

%% compute posterior estimates
alpha1 = (gamma .^ 2 ./ (sig1 .^ 2 + gamma .^ 2));
alpha2 = (gamma .^ 2 ./ (sig2 .^ 2 + gamma .^ 2));
p1 = vel1 .* alpha1;
p2 = vel2 .* alpha2;

%% choices are maximum of those
T = p2 > p1;

%% save for later use
s3 = [v1 c1 v2 c2 nan(nTrials, 4) T nan(nTrials, 2)];
save('/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/NewCode/Data/s3.mat', ...
    's3', 'sig1', 'sig2');

% save information about distributions
contrasts = testContrasts;
sigmas = (1 ./ testContrasts) .* maxSigma;
% sigmas = sigs;
vels = testVels;
save(['/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/NewCode/Data/' ...
    'gaussianDists.mat'], 'sigmas', 'vels', 'gamma', 'contrasts', 'sig1', 'sig2');