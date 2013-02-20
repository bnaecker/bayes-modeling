function ds = getStatistics(ds)
%
% function ds = getStatistics(ds)
%
% Computes statistics of bootstrapped datasets
%
% INPUT:    data - the main data structure
%
% OUTPUT:   data - same
%
% (c) bnaecker@stanford.edu 24 Mar 2012

%% Notify
fprintf('\nComputing statistics of bootstrapped fits...');

%% Compute statistics
for subject = 1:ds.info.nSubjects
    ds.info.currentSubject = subject;
    
    % Variance of likelihood
    ds.stats(subject).hVar = var(ds.params(subject).hHat);
    
    % Variance of prior
    ds.stats(subject).priorVar = var(ds.params(subject).prior);
    
    % Calculate +/- 2 SDs from mean
    ds.stats(subject).hLB = ds.params(subject).hHat(1, :) - ds.info.nSEs .* ...
        sqrt(ds.stats(subject).hVar);
    ds.stats(subject).hUB = ds.params(subject).hHat(1, :) + ds.info.nSEs .* ...
        sqrt(ds.stats(subject).hVar);
    
    ds.stats(subject).priorLB = max(1e-12, ds.params(subject).prior(1, :) - ...
        ds.info.nSEs .* sqrt(ds.stats(subject).priorVar));
    ds.stats(subject).priorUB = ds.params(subject).prior(1, :) + ...
        ds.info.nSEs .* sqrt(ds.stats(subject).priorVar);
    
    % If not fitting g(v) (speed-dependence of likelihood width), divide
    % h(c) by assumed g of each subject (pulled from figure 4 in Stocker &
    % Simoncelli)
    if ~ds.flags.fitG
        ds.stats(subject).likeLB = ds.params(subject).likeWidth(1, :) - ...
            ds.info.nSEs .* sqrt(ds.stats(subject).hVar);
        ds.stats(subject).likeUB = ds.params(subject).likeWidth(1, :) + ...
            ds.info.nSEs .* sqrt(ds.stats(subject).hVar);
    end
end

%% Organize
if ds.flags.fitLikeFun
    str = 'done.';
else
    str = 'done.\n';
end
fprintf(str);
ds = orderDataStruct(ds);