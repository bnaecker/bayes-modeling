function ds = getFitStatistics(ds)
%
% FUNCTION ds = getFitStatistics(ds)
%
% Computes statistics over the resampled data sets
%
% (c) bnaecker@stanford.edu 11 Nov 2012

%% notify
fprintf('\nComputing statistics of bootstrapped fits ... ');

%% compute statistics
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    % variance of likelihood 
    ds.stats(si).likeVar = var(ds.params(si).likeHat, 0, 2);
    
    % variance of prior
    ds.stats(si).priorVar = var(ds.params(si).prior, 0, 2);
    
    % calculate +/- requested standard errors
    ds.stats(si).likeLB = ds.params(si).likeHat(:, 1) - ...
        ds.info.nSEs .* sqrt(ds.stats(si).likeVar);
    ds.stats(si).likeUB = ds.params(si).likeHat(:, 1) + ...
        ds.info.nSEs .* sqrt(ds.stats(si).likeVar);
    
    ds.stats(si).priorLB = max(1e-12, ds.params(si).prior(:, 1) - ...
        ds.info.nSEs .* sqrt(ds.stats(si).priorVar));
    ds.stats(si).priorUB = ds.params(si).prior(:, 1) + ...
        ds.info.nSEs .* sqrt(ds.stats(si).priorVar);
    
    % get the bounds to the likelihood widthds, depending on whether or not
    % we've fit the speed-dependence of the likelihods
    if ds.flags.fitLikeSpeed
        % optional functionality
    else
        ds.stats(si).likeWidthLB = ds.params(si).likeWidth(:, 1) - ...
            ds.info.nSEs .* sqrt(ds.stats(si).likeVar);
        ds.stats(si).likeWidthUB = ds.params(si).likeWidth(:, 1) + ...
            ds.info.nSEs .* sqrt(ds.stats(si).likeVar);
    end
end

%% organize
fprintf('done.\n');
ds = orderDataStruct(ds);
