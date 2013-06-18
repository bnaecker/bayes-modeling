function sim = simulateDiscriminationData(varargin)
%
% function sim = simulateDiscriminationData(varargin)
%
% Helper function that simulates a Bayesian ideal observer using the
% requested prior. The output is 2AFC data which can then be fitted using
% the included functions. This is largely used to validate the fitting
% procdures.
%
% (c) bnaecker@stanford.edu 13 Nov 2012

%% check varargin
% prior type and definition
if any(strcmpi('priortype', varargin))
    sim.priorType = varargin{find(strcmpi('priortype', varargin)) + 1};
    if any(strcmpi('gamma', varargin))
        sim.gamma = varargin{find(strcmpi('gamma', varargin)) + 1};
    else
        sim.gamma = sqrt(10);
    end
    if any(strcmpi('slopes', varargin))
        sim.priorSlopes = varargin{find(strcmpi('slopes', varargin)) + 1};
    else
        sim.priorSlopes = [6.5 8.5 9 7.9 2 0.05]';
    end
else
    sim.priorType = 'loglinear';
    if any(strcmpi('slopes', varargin))
        sim.priorSlopes = varargin{find(strcmpi('slopes', varargin)) + 1};
    else
        sim.priorSlopes = [6.5 8.5 9 7.9 2 0.05]';
    end
end
% check priorType is something we can deal with
assert(any(strcmpi(sim.priorType, {'gaussian', 'loglinear'})), ...
    'simulateDiscriminationData:unknownPriorType', ...
    'Supported prior types are "gaussian" and "loglinear"');

% velocities
if any(strcmpi('vels', varargin))
    sim.vels = varargin{find(strcmpi('vels', varargin)) + 1};
else
    sim.vels = [0.5 1 2 4 8 12]';
end

% contrasts
if any(strcmpi('contrasts', varargin))
    sim.contrasts = varargin{find(strcmpi('contrasts', varargin)) + 1};
else
    sim.contrasts = [0.05 0.075 0.1 0.2 0.4 0.5 0.8]';
    c1IndVals = [1 7];
end

%% notify
fprintf('\nSimulating ideal observer with %s prior ... ', sim.priorType);

%% setup the trials
% number of trials
sim.nTrials = 1e4;

% contrasts
sim.c1Inds = randi(c1IndVals, sim.nTrials, 1);
sim.c2Inds = randi(length(sim.contrasts), sim.nTrials, 1);
sim.c1 = sim.contrasts(sim.c1Inds);
sim.c2 = sim.contrasts(sim.c2Inds);

% "reference" velocities
v1Inds = randi(length(sim.vels), sim.nTrials, 1);
sim.v1 = sim.vels(v1Inds);

% approximate an adaptive staircase procedure for test velocities
minF = 0.8; maxF = 1.2;
nF = 10;
f = linspace(minF, maxF, nF);
testV = sim.vels * f;
sigma = 1;
weight = normpdf(linspace(-1, 1, nF), 0, sigma);
weight = weight ./ sum(weight);
sim.v2 = zeros(sim.nTrials, 1);
for ti = 1:sim.nTrials
    vsample = testV(v1Inds(ti), :);
    sim.v2(ti) = randsample(vsample, 1, true, weight);
end

%% compute the values drawn from the likelihood on each trial
% functional form for likelihood width
sim.q = 1.4;
sim.c50 = .16;
sim.rMin = 0.45;
sim.rMax = 3;
sim.likeWidth = 1 ./ sqrt(sim.rMax .* (sim.contrasts .^ sim.q ./ ...
    (sim.contrasts .^ sim.q + sim.contrasts .^ sim.c50)) + sim.rMin);

% add gaussian noise to presented velocities (draw from likelihood
% distributions)
sim.vel1 = max(sim.v1 + randn(sim.nTrials, 1) .* sim.likeWidth(sim.c1Inds), 0);
sim.vel2 = max(sim.v2 + randn(sim.nTrials, 1) .* sim.likeWidth(sim.c2Inds), 0);

%% make prior distribution
sim.dv = 0.01;
sim.vAx = min(sim.vels):sim.dv:(max(sim.vels) + 0.5 * max(sim.vels));
if strcmp(sim.priorType, 'gaussian')
    
    % gaussian prior
    sim.prior = normpdf(sim.vAx, 0, sim.gamma);
    sim.prior = sim.prior ./ (sum(sim.prior) * sim.dv);
    
    % look up the value of the prior at the tested velocities
    sim.prior1 = interp1(sim.vAx, sim.prior, sim.vel1, 'linear', 'extrap');
    sim.prior2 = interp1(sim.vAx, sim.prior, sim.vel2, 'linear', 'extrap');
elseif strcmp(sim.priorType, 'loglinear')
    
    % loglinear prior, determined by its slope
    sim.allPriorSlopes = interp1(sim.vels, sim.priorSlopes, sim.vAx, ...
        'linear', 'extrap');
    sim.prior = exp(-cumsum(sim.allPriorSlopes) * sim.dv);
    sim.prior = sim.prior ./ (sum(sim.prior) * sim.dv);
    
    % lookup the prior slope at the tested velocities
    sim.prior1 = interp1(sim.vels, sim.priorSlopes, sim.v1, 'linear', 'extrap');
    sim.prior2 = interp1(sim.vels, sim.priorSlopes, sim.v2, 'linear', 'extrap');
end

%% compute posterior estimates
if strcmp(sim.priorType, 'gaussian')
    sim.p1 = sim.vel1 .* ...
        (sim.gamma .^ 2 ./ (sim.likeWidth(sim.c1Inds) .^ 2 + sim.gamma .^ 2));
    sim.p2 = sim.vel2 .* ...
        (sim.gamma .^ 2 ./ (sim.likeWidth(sim.c2Inds) .^ 2 + sim.gamma .^ 2));
elseif strcmp(sim.priorType, 'loglinear')
    sim.p1 = sim.vel1 - sim.prior1 .* sim.likeWidth(sim.c1Inds) .^ 2;
    sim.p2 = sim.vel2 - sim.prior2 .* sim.likeWidth(sim.c2Inds) .^ 2;
end

%% make choices
sim.T = sim.p2 > sim.p1;

% make lapse rate
lapse = .01;
swp = rand(sim.nTrials, 1) < lapse;
sim.T(swp) = ~sim.T(swp);

%% save data
% notify
fprintf('done.\nSaving ... ');
m = [sim.v1 sim.c1 sim.v2 sim.c2 nan(sim.nTrials, 4) sim.T nan(sim.nTrials, 2)];
save(['../../data/simData_' sim.priorType '.mat'], 'm');
fprintf('done.\n');
