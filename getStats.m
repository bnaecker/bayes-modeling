function var = getStats(var)
%
% getStats calculates the statistics of bootstrapped data sets
%
% (c) benjamin naecker UT Austin 26 May 2011 benjamin.naecker@gmail.com

% For each subject
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    % Check that these are all within var.info.SDThresh. If not, refit.
    good = 0;
%     while ~good
%         [var, good] = checkFits(var);
%     end
    
    % Get variance of likelihood
    var.stats(subject).hStd = std(var.params(subject).hHat);
    
    % Get variance of prior
    var.stats(subject).prStd = std(var.params(subject).prior);
    
    % Use these to calculate +/- 2 SEs from the mean
    var.params(subject).hLB = var.params(subject).hHat(1, :) - var.info.nSEs .* ...
        var.stats(subject).hStd;
    var.params(subject).hUB = var.params(subject).hHat(1, :) + var.info.nSEs .* ...
        var.stats(subject).hStd;
    
    var.params(subject).prLB = max(1e-12, var.params(subject).prior(1, :) - var.info.nSEs .* ...
        var.stats(subject).prStd);
    var.params(subject).prUB = var.params(subject).prior(1, :) + var.info.nSEs .* ...
        var.stats(subject).prStd;
    
    % If we're not fitting g(v), divide h(C) by the assumed g of each
    % subject
    if ~var.flags.fitG
        var.params(subject).likeLB = var.params(subject).likeWidth(1, :) - ...
            var.info.nSEs .* var.stats(subject).hStd;
        var.params(subject).likeUB = var.params(subject).likeWidth(1, :) + ...
            var.info.nSEs .* var.stats(subject).hStd;
    end
end

%% Organize
var = organizeVar(var);