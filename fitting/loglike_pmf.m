function L = loglike_pmf(prs, fitType, X, choice)
%
% helper function computes the log-likelihood of the weibull fit to the
% data given the parameters. May be updated to allow for fitting a
% cumulative normal to the data instead of a weibull.
%
% (c) bnaecker@stanford.edu 14 Nov 2012

%% Compute the CDF with the current parameters
% fit = cdf(fitType, X, prs(1), prs(2));
fit = 1 - exp(-(X ./ prs(1)) .^ prs(2));

%% Bound values away from 0 and 1
EPS = 1e-12;
fit = max(min(fit, 1 - EPS), EPS);

%% Compute log-likelihood
L = sum(choice .* log(fit) + (1 - choice) .* (1 - fit));