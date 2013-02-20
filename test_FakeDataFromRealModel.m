%% Script: Runs the speed disrcimiantion model with the requested prior
% on real data, uses the fitted parameters to generate fake data, and then
% fits PMFS to that data for comparison with the model generated PMFs.

% (c) bnaecker@stanford.edu 18 Jun 2012

%% Run model on requested prior type
fprintf('\nFitting model to real data...\n')
priorType = 'loglinear';
nSubjects = 1;
ds = runSpeedDiscriminationModel('nSubjects', nSubjects, 'priorType', priorType);

%% Convert log-transformed parameters back to linear coordinates
fprintf('\nSimulating data...');
refVels = ds.velTrans.iTransFun(ds.data.refVels);

%% Generate fake data with these parameters
% Contrasts
nTrials = 1e4;
c1inds = randi(length(ds.data.refContrasts), nTrials, 1);
c2inds = randi(length(ds.data.testContrasts), nTrials, 1);
c1 = ds.data.refContrasts(c1inds);
c2 = ds.data.testContrasts(c2inds);

likeWidths = ds(1).params.likeWidth(:, 1);
lowLikeWidths = [likeWidths(2) likeWidths(6)];
sig1 = likeWidths(c1inds);
sig2 = likeWidths(c2inds);

% Speeds (test vels are randomly sampled from interval around ref vel)
v1inds = randi(length(refVels), nTrials, 1);
v1 = refVels(v1inds);

%% Make test velocities
minFrac = .5; maxFrac = 1.5;
nFracs = 5;
f = linspace(minFrac, maxFrac, nFracs);
f = ones(length(refVels), 1) * f;
testV = refVels * ones(1, nFracs);
testV = testV .* f;

%% Sample these velocities appropriately
sigma = 1;
weight = normpdf(linspace(-.5, .5, nFracs), 0, sigma);
weight = weight ./ sum(weight);
v2 = zeros(size(v1));
for ti = 1:nTrials
    vs = testV(v1inds(ti), :);
    v2(ti) = randsample(vs, 1, true, weight);
end

%% Add Gaussian noise to the samples from the likelihood widths
vel1 = v1 + randn(nTrials, 1) .* sig1;
vel2 = v2 + randn(nTrials, 1) .* sig2;

%% Compute correct prior description from the fitted values
priorSlopes = ds.params.prior(:, 1);
prior1 = priorSlopes(v1inds);
prior2 = interp1(refVels, priorSlopes, vel2, 'linear', 'extrap');

%% Compute posterior estimates
alpha1 = prior1 .* sig1 .^ 2;
alpha2 = prior2 .* sig2 .^ 2;
p1 = vel1 - alpha1;
p2 = vel2 - alpha2;

%% Simulate choices. MAP estimates
T = p2 > p1;
fprintf('\nSaving data...\n');
s4 = [v1 c1 v2 c2 nan(nTrials, 4) T nan(nTrials, 2)];
save('/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/Code/Data/s4.mat', ...
    's4', 'sig1', 'sig2');

% % save information about distributions, for plotting/comparisons
% contrasts = testContrasts;
% simgas = likeWidths;
% vels = testVels;
% save(['/Users/bnaecker/FileCabinet/Projects/SpeedDiscrimination/NewCode/Data/' ...
%     'loglinearDists.mat'], 'sigmas', 'vels', 'priorSlopes', 'contrasts', 'sig1', 'sig2');

%% 
fprintf('Fitting new model to simulated data...\n');
newData = runSpeedDiscriminationModel('nSubjects', 0, 'priorType', priorType);
% newData = plotPsychometricFunctions(newData);