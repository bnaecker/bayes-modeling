function L = logli_var(prs, var)
%
% logli_var takes the input struct 'var' and returns the log-likelihood of
% the current parameters given the data
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

% Indices of current subject and bootstrap
subject = var.info.currentSubject;
boot = var.info.currentBoot;

if boot > 1
    stophere = 1;
end

% To bound LL away from zero later
EPS = 1e-12;

% Assign current parameters guessed from fmincon
sigVals = prs(1:var.info.nUniqueContrasts);
a = prs(var.info.nUniqueContrasts+1:end);

% Arrange slopes of prior at each given speed
refSlope = interp1(var.data(subject).refVels, a,...
                   var.data(subject).refV, 'linear');
testSlope = interp1(var.data(subject).refVels, a,...
                    var.data(subject).testV, 'linear', 'extrap');
                
% Arrange likelihood widths h(c)
refSig = sigVals(var.data(subject).refCInds);
testSig = sigVals(var.data(subject).testCInds);

% Estimate distribution
Z = (var.data(subject).testV - var.data(subject).refV + refSlope .* refSig .^2 ...
    - testSlope .* testSig .^ 2) ./...
    sqrt(refSig .^ 2 + testSig .^ 2);

% Normcdf of estimate distribution, bounded away from zero
rho = min(max(normcdf(Z), EPS), 1 - EPS);

% Calculate log likelihood
L = sum(var.data(subject).T(:, boot) .* log(rho) + (1 - var.data(subject).T(:, boot)) .* log(1 - rho));