function err = fitCDF(prs, fitType, X, choice)
%
% helper function fits either 'wbl' or 'normal' cdf to the psychometric
% data from speed discrimination experiments. Maximum-likelihood fitting.
%
% (c) benjamin naecker UT Austin 31 May 2011 benjamin.naecker@gmail.com

%% Compute the CDF with the current parameters
% fit = cdf(fitType, X, prs(1), prs(2));
fit = 1 - exp(-(X ./ prs(1)) .^ prs(2));

%% Bound values away from 0 and 1
EPS = 1e-12;
fit = max(min(fit, 1 - EPS), EPS);

%% Compute log-likelihood
err = -sum(choice .* log(fit) + (1 - choice) .* (1 - fit));