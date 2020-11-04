function ds = computePsychometricFunctions(ds)
%
% FUNCTION ds = computePsychometricFunctions(ds)
%
% Computes the psychometric functions for the Bayesian ideal observer model,
% both those fitted to data and those actually predicted from the model
% itself.
%
% (c) bnaecker@stanford.edu 14 Nov 2012

% loop over subjects
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    % preallocate a cell array to hold the proportion of times the test
    % grating was seen faster, the fitted and model-predicted pmfs, and the
    % domain over which these are defined
    ds.pmfs(si).proportionTestFaster = cell(ds.info.nUniqueRefContrasts, ...
        ds.info.nUniqueRefVels, ds.info.nUniqueContrasts);
    ds.pmfs(si).pmfPrs = ds.pmfs(si).proportionTestFaster;
    ds.pmfs(si).dataPmfs = ds.pmfs(si).proportionTestFaster;
    ds.pmfs(si).modelPmfs = ds.pmfs(si).proportionTestFaster;
    ds.pmfs(si).pmfDomain = ds.pmfs(si).proportionTestFaster;
    
    % loop over each unique triplet of conditions, starting with reference contrast
	% notify 
	fprintf('  computing pmfs ');
    for refCNum = 1:ds.info.nUniqueRefContrasts
        % loop over reference velocity
        for refVNum = 1:ds.info.nUniqueRefVels
            % loop over test contrast
            for testCNum = 1:ds.info.nUniqueContrasts
				% notify
				if (mod((refCNum - 1) * ds.info.nUniqueRefVels * ds.info.nUniqueContrasts + ...
					(refVNum - 1) * ds.info.nUniqueContrasts + testCNum, ...
					floor(ds.info.nUniqueRefContrasts * ds.info.nUniqueRefVels * ...
					ds.info.nUniqueContrasts / 10)) == 0)
					fprintf('.');
				end
				%fprintf('  computing pmf for condition %d of %d\n', ...
					%(refCNum - 1) * ds.info.nUniqueRefVels * ds.info.nUniqueContrasts + ...
					%(refVNum - 1) * ds.info.nUniqueContrasts + ...
					%testCNum, ...
					%ds.info.nUniqueRefContrasts * ds.info.nUniqueRefVels * ...
					%ds.info.nUniqueContrasts);
                
                % indices into this triplet condition
                refCInds = ds.data(si).refC == ...
                    ds.data(si).refContrasts(refCNum);
                refVInds = ds.data(si).refV == ...
                    ds.data(si).refVels(refVNum);
                testCInds = ds.data(si).testC == ...
                    ds.data(si).testContrasts(testCNum);
                
                % get the test and reference velocities for this triplet condition
                testVels = unique(ds.data(si).testV(...
                        refCInds & refVInds & testCInds));
                refVels = ds.data(si).refVels; 
                
                % check that there is data for this condition
                if ~isempty(testVels)
                    
                    % preallocate the array to hold the proportion of times
                    % the test grating was reported as moving faster
                    ds.pmfs(si).proportionTestFaster{refCNum, refVNum, testCNum} = ...
                        nan(1, length(testVels));
                    
                    % grab the speed and choice data explicitly, so that we 
                    % can do maximum-likelihood fitting of the pmf
					speed = ds.data(si).testV(refCInds & refVInds & testCInds);
                    choice = ds.data(si).T(refCInds & refVInds & testCInds);
                    
                    % loop over each test velocity presented at this
                    % triplet condition, computing the proportion of times
                    % the test grating was reported as moving faster
                    for testVNum = 1:length(testVels)
                        
                        % indices to the choices made on trials with this
                        % particular test velocity
                        choiceInds = refCInds & refVInds & testCInds & ...
                            ds.data(si).testV == testVels(testVNum);
                        
                        % calculate proportion of times in which the test
                        % grating was reported faster, i.e., choice == 1
                        ds.pmfs(si).proportionTestFaster{refCNum, ...
                            refVNum, testCNum}(1, testVNum) = ...
                            sum(ds.data(si).T(choiceInds, 1)) ./ ...
                            length(ds.data(si).T(choiceInds, 1));
                    end
                    
                    % fit a cumulative weibull function to the  choice data
                    if ds.flags.useCFit
                        
                        % define the function to be fit
                        weibullFun = '1 - exp(-(x ./ lambda) .^ k)';
                        
                        % set a start point for the minimization
                        % [scale; shape]
                        startPoint = [refVels(refVNum); ...
                            1 / ds.data(si).testContrasts(testCNum)];
                        
                        % fit options
                        fitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                            'StartPoint', startPoint, 'Lower', [0 0]);
                        
                        % create the fit object
                        fitObj = fittype(weibullFun, 'Independent', 'x', ...
                            'Coefficients', {'lambda', 'k'}, 'Options', fitOpt);
                        
                        % use the curve fitting toolbox to create the fit
                        try
                        theFit = fit(testVels, ...
                            ds.pmfs(si).proportionTestFaster{refCNum, refVNum, ...
                            testCNum}', fitObj);
                        catch me
                            h = 1;
                        end
                        
                        % save the parameters
                        ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum} = ...
                            [theFit.k theFit.lambda];
                        
                        % generate the domain of the pmf
                        ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum} = ...
                            linspace(0.5 * refVels(refVNum) , ...
                            2 * refVels(refVNum), 1000);
                        
                        % use the fit object to create the pmf
                        ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum} = ...
                            theFit(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum});
                    else
                        % no curve fitting toolbox, so we'll fit the pmf by
                        % maximum-likelihood
                        
                        % initial parameter values, with bounds
                        prs0 = [refVels(refVNum); ...
                            ds.data(si).testContrasts(testCNum)];
                        LB = [1e-15; 1e-15];
                        UB = [100; 100];
                        
                        % define the objective function, the negative
                        % log-likelihood of the the fit
                        fitFun = @(prs) -loglike_pmf(prs, ds.flags.pmfFitType, ...
                            speed, choice);
                        
                        % minimize the negative log-likelihood of the fit
                        % using fmincon
                        ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum} = ...
                            fmincon(fitFun, prs0, [], [], [], [], LB, UB, ...
                            [], ds.llOpt.options);
                        
                        % define the domain of the pmf
                        ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum} = ...
                            linspace(0.5 * refVels(refVNum), ...
                            2 * refVels(refVNum), 1000);
                        
                        % use the fitted parameters to compute the pmf
                        ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum} = ...
                            1 - exp(-(ds.pmfs(si).pmfDomain{refCNum, refVNum, ...
                            testCNum} ./ ds.pmfs(si).pmfPrs{refCNum, refVNum, ...
                            testCNum}(1)) .^ ds.pmfs(si).pmfPrs{refCNum, refVNum, ...
                            testCNum}(2));
                    end
                    
                    % compute the model-predicted pmfs
                    if strcmp(ds.flags.priorType, 'gaussian')

                        % grab width of the prior (same for all velocities)
                        % and the width of the likelihood for the current
                        % condition, first for the reference stimulus
                        x1 = refVels(refVNum);
                        sig1 = ds.params(si).likeWidth(...
                            ds.data(si).refContrasts(refCNum) == ...
                            ds.data(si).testContrasts, 1);
                        sig2 = ds.params(si).likeWidth(testCNum, 1);
                        
                        % now for the test stimulus
                        x2 = ds.pmfs(si).pmfsDomain{refCNum, refVNum, testCNum};
                        
                        % compute alpha1 and alpha2, see Naecker & Pillow
                        % for full description. opposite from the sigmas
                        a1 = ds.params(si).gamma(1) .^ 2 + sig2 .^ 2;
                        a2 = ds.params(si).gamma(1) .^ 2 + sig1 .^ 2;
                        
                        % compute the predicted psychometric function.
                        % again, see Naecker & Pillow for full description
                        ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum} = ...
                            sort(normcdf((a2 .* x2 - a1 .* x1) ./ ...
                            sqrt(a1 .^ 2 .* sig1 .^2 + a2 .^ 2 .* sig2 .^ 2)));
                    elseif strcmp(ds.flags.priorType, 'loglinear')
                        
                        % grab the slope of the log-prior and the
                        % likelihood width for the current condition, first
                        % of the reference stimulus
                        %x1 = refVels(refVNum);
                        %a1 = ds.params(si).slopeHat(refVNum, 1);
                        %sig1 = ds.params(si).likeWidth(...
                            %ds.data(si).refContrasts(refCNum) == ...
                            %ds.data(si).testContrasts, 1);
                        %sig2 = ds.params(si).likeWidth(testCNum, 1);
                        
                        %% then for the test stimulus
                        %x2 = ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum};
                        %a2 = interp1(refVels, ...
                            %ds.params(si).slopeHat(:, 1), x1, 'linear', 'extrap');
                        
                        %% compute the predicted psychometric function
                        %ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum} = ...
                            %sort(normcdf(((x2 - x1 + a1 .* sig1 .^ 2 - ...
                            %a2 .* sig2 .^ 2) ./ sqrt(sig1 .^ 2 + sig2 .^ 2))));
						
						% simulate data for each point on the pmf axis
						for vi = 1:length(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum})
							% make the reference measurements (draws from encoding distribution)
							refMeas = max(refVels(refVNum) + ...
								ds.params(si).likeWidth(refCNum, 1) .* ...
								randn(1, ds.flags.nPmfSamples), 0);

							% make the test measurements (draws from encoding distributions)
							testMeas = max(...
								ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}(vi) + ...
								ds.params(si).likeWidth(refCNum, 1) .* ...
								randn(1, ds.flags.nPmfSamples), 0);

							% make the reference likelihood function
							likeSupport = (ones(ds.flags.nPmfSamples, 1) * ...
								ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum})';
							refLikeMeans = ones(length(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}), 1) * ...
								refMeas;
							refLike = normpdf(likeSupport, refLikeMeans, ...
								ds.params(si).likeWidth(refCNum, 1));

							% make the test likelihood function
							testLikeMeans = ones(1, length(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum})) * ...
								testMeas;
							testLike = normpdf(likeSupport, testLikeMeans, ...
								ds.params(si).likeWidth(testCNum, 1));

							% get the the prior distributions
							refPrior = interp1(ds.params(si).interpAx, ...
								ds.params(si).interpPrior(:, 1), ...
								refMeas, 'linear', 'extrap');
							refPrior = refPrior ./ (sum(refPrior) * ds.params(si).priorDx);
							testPrior = interp1(ds.params(si).interpAx, ...
								ds.params(si).interpPrior(:, 1), ...
								testMeas, 'linear', 'extrap');
							testPrior = testPrior ./ (sum(testPrior) * ds.params(si).priorDx);

							% compute the posterior distributions
							refPost = refMeas .* refPrior;
							refPost = refPost ./ (sum(refPost) * ds.params(si).priorDx);
							testPost = testMeas .* testPrior;
							testPost = testPost ./ (sum(testPost) * ds.params(si).priorDx);

							% compare the posterior estimates, i.e., fraction of times the
							% test was 'seen' as faster
							ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}(vi) = ...
								sum(testPost > refPost) / ds.flags.nPmfSamples;
						end
						h = 1;
                    else
                        % optional functionality
                    end
                end
            end
        end
    end
	fprintf('done.\n');
end
