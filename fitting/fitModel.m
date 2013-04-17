function ds = fitModel(ds)
%
% FUNCTION ds = fitModel(ds)
%
% This function fits the model assuming the requested prior type to the
% psychophysical discrimination data in the data structure ds.
%
% (C) bnaecker@stanford.edu 11 Nov 2012

%% fit model for each subject
for si = 1:ds.info.nSubjects
    % notify
    fprintf('\nFitting %s model for subject %d\n\n', ds.flags.priorType, si);
    
    % setup indices
    ds.info.currentSubject = si;
    ds.info.currentBoot = 1;
    
    % define all the unique triplet conditions (refV, refC, testC). use
    % this to index into choices during the resampling procedure
    conditions = cartprod(ds.data(si).refVels, ...
        ds.data(si).refContrasts, ds.data(si).testContrasts);
    
    % preallocate arrays to hold output of minimization procedure, for each
    % resampled data set
    ds.params(si).likeHat = nan(ds.info.nUniqueContrasts, ds.info.nBoots);
    ds.params(si).likeWidth = nan(ds.info.nUniqueContrasts, ds.info.nBoots);
	switch ds.flags.priorType
    	case 'loglinear'
        	ds.params(si).slopeHat = nan(ds.info.nUniqueRefVels, ds.info.nBoots);
        	ds.llOpt(si).prsHat = nan(ds.info.nUniqueContrasts + ds.info.nUniqueRefVels, ...
        	ds.info.nBoots);
    	case 'gaussian'
        	ds.params(si).gamma = nan(ds.info.nBoots, 1);
        	ds.llOpt(si).prsHat = nan(ds.info.nUniqueContrasts + 1, ds.info.nBoots);
		case 'mixture'
			ds.params(si).mixMean = nan(ds.mixInfo.nComponents, ds.info.nBoots);
			ds.params(si).mixVar = nan(ds.mixInfo.nComponents, ds.info.nBoots);
			ds.llOpt(si).prsHat = nan(ds.info.nUniqueContrasts + ...
				2 .* ds.mixInfo.nComponents, ds.info.nBoots);
    end
    ds.llOpt(si).fval = nan(ds.info.nBoots, 1);
    ds.llOpt(si).eflag = nan(ds.info.nBoots, 1);
    
    % outermost fitting loop, over bootstraps
    while ds.info.currentBoot <= ds.info.nBoots
        
        % resample appropriately, i.e., choices among unique sets of conditions
        if ds.info.currentBoot > 1
            for ci = 1:size(conditions, 1)
                % get the current condition and its indices
                currentCondition = conditions(ci, :);
                conditionInds = ...
                    ds.data(si).refV == currentCondition(1) & ...
                    ds.data(si).refC == currentCondition(2) & ...
                    ds.data(si).testC == currentCondition(3);
                
                % make sure there are trials of this condition
                if sum(conditionInds) > 0
                    % find unique test velocities presented for this
                    % triplet
                    uniqueVels = unique(ds.data(si).testV(conditionInds));
                    
                    % get the subject's choices for this triplet
                    d = ds.data(si).T(conditionInds, 1);
                    
                    % preallocated resample array
                    reSample = nan(size(d));
                    
                    % loop over unique test velocities in this triplet
                    for vi = 1:length(uniqueVels)
                        
                        % find indices of this velocity
                        vInds = uniqueVels(vi) == ...
                            ds.data(si).testV(conditionInds);
                        
                        % deal with case where there is only one of this
                        % velocity
                        if sum(vInds) == 1
                            reSample(vInds) = d(vInds);
                        else
                            % if there are multiple, randomly resample data
                            % with replacement
                            reSample(vInds) = ...
                                randsample(d(vInds), sum(vInds), true)';
                        end
                    end
                    
                    % allocate the resampled data to the full choice array
                    ds.data(si).T(conditionInds, ds.info.currentBoot) = ...
                        reSample;
                end
            end
        end
        
        % notify which bootstrap is being fit currently
        fprintf('Boot %d of %d ... ', ds.info.currentBoot, ds.info.nBoots);
        
        % redefine handle to log-likelihood function. this must be done
        % because we have updated the data structure, and so must redefine
        % the anonymous function handle as well.
		switch ds.flags.priorType
        	case 'loglinear'
            	ds.llOpt(si).llFun = ...
                	@(prs) (loglike_loglinear(prs, ds));
        	case 'gaussian'
            	ds.llOpt(si).llFun = ...
                	@(prs) (loglike_gaussian(prs, ds));
			case 'mixture'
				ds.llOpt(si).llFun = ...
					@(prs) (loglike_mixture(prs, ds));
        	otherwise
            	% optional functionality
        end
        
        % actually run the minimization
        [ds.llOpt(si).prsHat(:, ds.info.currentBoot) ...
         ds.llOpt(si).fval(ds.info.currentBoot) ...
         ds.llOpt(si).eflag(ds.info.currentBoot)] = ...
            fmincon(ds.llOpt(si).llFun, ds.llOpt(si).prs0,...
            [], [], [], [], ds.llOpt(si).LB, ds.llOpt(si).UB, ...
            [], ds.llOpt(si).options);
        
        % save data, if the fitting converged. if the fitting does not
        % converge, the user must decide what should be done. for example,
        % one might rerun the optimization from a different starting point.
        switch ds.llOpt(si).eflag(ds.info.currentBoot)
            case 0
                % number of iterations exceeded either MaxIter or
                % MaxFunEvals
            case -2
                % no feasible point found
            otherwise 
                % notify
                fprintf('converged @ %.2f\n', ds.llOpt(si).fval(ds.info.currentBoot));
                
                % save the parameters
                ds.params(si).likeHat(:, ds.info.currentBoot) = ...
                    ds.llOpt(si).prsHat(1:ds.info.nUniqueContrasts, ...
                    ds.info.currentBoot);
                switch ds.flags.priorType
					case 'loglinear'
                    	ds.params(si).slopeHat(:, ds.info.currentBoot) = ...
                        	ds.llOpt(si).prsHat(ds.info.nUniqueContrasts + 1 : end, ...
                        	ds.info.currentBoot);
                	case 'gaussian'
                    	ds.params(si).gamma(ds.info.currentBoot) = ...
                        	ds.llOpt(si).prsHat(end, ds.info.currentBoot);
					case 'mixture'
						ds.params(si).mixMean(:, ds.info.currentBoot) = ...
							ds.llOpt(si).prsHat(ds.info.nUniqueContrasts + 1 : ...
							ds.info.nUniqueContrasts + ds.mixInfo.nComponents, ...
							ds.info.currentBoot);
						ds.params(si).mixVar(:, ds.info.currentBoot) = ...
							ds.llOpt(si).prsHat(ds.info.nUniqueContrasts + ...	
							ds.mixInfo.nComponents + 1 : end, ds.info.currentBoot);
					otherwise
                    	% optional functionality
                end
                
                % must include the speed-dependence of the likelihood
                % width, if it is not being fit
                if ds.flags.fitLikeSpeed
                    % optional functionality
                else
                    ds.params(si).likeWidth(:, ds.info.currentBoot) = ...
                        ds.params(si).likeHat(:, ds.info.currentBoot) ./ ...
                        ds.data(si).G;
                end
                
                % reconstruct the prior from the fitted prior slopes or
                % variances
                ds = reconstructPrior(ds);
                
                % increment boot counter
                ds.info.currentBoot = ds.info.currentBoot + 1;
        end
    end
end

%% order
ds = orderDataStruct(ds);
