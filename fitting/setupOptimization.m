function ds = setupOptimization(ds)
%
% FUNCTION ds = setupOptimization(ds)
%
% Setup optimization of model log-likelihood function.
%
% (c) bnaecker@stanford.edu 11 Nov 2012

fprintf('Setting up optimization parameters ... ');
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    %% setup minimization of log-likelihood of model
    % function handle to log-likelihood
    if ds.flags.fitLikeSpeed
        % optional functionality
	else
		switch ds.flags.priorType
			case 'loglinear'
        		ds.llOpt(si).llFun = @(prs) (loglike_loglinear(prs, ds));
    		case 'gaussian'
        		ds.llOpt(si).llFun = @(prs) (loglike_gaussian(prs, ds));
			case 'mixture'
				ds.llOpt(si).llFun = @(prs) (loglike_mixutre(prs, ds));
    		otherwise
        		error('runSpeedDiscriminationModel:setupOptimization:unknownPrior', ...
            		'The supported prior types are "loglinear" or "gaussian"');
		end
    end
    
    %% options to fmincon
   	ds.llOpt(si).options = ...
    	optimset('Algorithm', 'interior-point', 'Display', 'off', ...
        'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');
    
    %% initial values for likelihood widths and prior parameters
	switch ds.flags.priorType
		case 'loglinear'
			% guesses and bounds to likelihood widths
        	likeGuess =  1 .* ones(ds.info.nUniqueContrasts, 1);
        	likeLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        	likeUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        
        	% guesses and bounds to prior slopes
        	priorGuess = 2 .* ones(ds.info.nUniqueRefVels, 1);
        	priorLB = -100 .* ones(ds.info.nUniqueRefVels, 1);
        	priorUB = 100 .* ones(ds.info.nUniqueRefVels, 1);

			% mixture weights are empty
			mixWeightsGuess = [];
			mixWeightsLB = [];
			mixWeightsUB = [];
        
    	case 'gaussian'
        	% guesses and bounds to likelihood widths
        	likeGuess = 1 .* ones(ds.info.nUniqueContrasts, 1);
        	likeLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        	likeUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        
        	% guesses and bounds to prior variance
        	priorGuess = 10;
        	priorLB = 1e-5;
        	priorUB = 1e5;

			% mixture weights are empty
			mixWeightsGuess = [];
			mixWeightsLB = [];
			mixWeightsUB = [];
        
    	case 'mixture'
			% guesses and bounds to likelihood widths
        	likeGuess = 1 .* ones(ds.info.nUniqueContrasts, 1);
        	likeLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        	likeUB = 100 .* ones(ds.info.nUniqueContrasts, 1);

			% guesses and bounds to mixture means
			priorGuess = [ds.mixInfo.mixMeans; ds.mixInfo.mixVars];
			priorLB = [-1e2 .* ones(ds.mixInfo.nComponents, 1); ...
					   1e-3 .* ones(ds.mixInfo.nComponents, 1)];
			priorUB = [1e2 .* ones(ds.mixInfo.nComponents, 1); ...
					   1e3 .* ones(ds.mixInfo.nComponents, 1)];

			% mixture weights
			mixWeightsGuess = ds.mixInfo.mixWeights;
			mixWeightsLB = zeros(ds.mixInfo.nComponents, 1);
			mixWeightsUB = ones(ds.mixInfo.nComponents, 1);
        
		otherwise
        	% optional functionality
    end
    
    % organize the guesses
    ds.llOpt(si).prs0 = [likeGuess; priorGuess; mixWeightsGuess];
    ds.llOpt(si).LB = [likeLB; priorLB; mixWeightsLB];
    ds.llOpt(si).UB = [likeUB; priorUB; mixWeightsUB];
    
    
    %% setup fitting of functional form to likelihood widths
    if ds.flags.fitLikeFun
        
        % initial values, including rates if requested
        if ds.flags.fitRates
            % [slope; semi-saturation point; minimum rate; maximum rate]
            ds.likeOpt(si).prs0 = [2; 0.2; 0.5; 5];
            ds.likeOpt(si).LB = [1e-15; 0.01; 0; 1];
            ds.likeOpt(si).UB = [100; 1; 20; 200];
        else
            % [slope; semi-saturation point]
            ds.likeOpt(si).prs0 = [2; 0.2];
            ds.likeOpt(si).LB = [1e-15; 1e-15];
            ds.likeOpt(si).UB = [100; 1];
            
            % min and max firing rates, if not fitted
            ds.likeOpt(si).R = [0.1; 100]; 
        end
        
        % use curve-fitting toolbox, if available
        v = ver('curvefit');
        ds.flags.useCFit = ~isempty(v);
        if ds.flags.useCFit;
            % fit options
            ds.likeOpt(si).fitOptions = ...
                fitoptions('Method', 'NonlinearLeastSquares', ...
                'Lower', ds.likeOpt(si).LB, 'Upper', ds.likeOpt(si).UB, ...
                'StartPoint', ds.likeOpt(si).prs0);
            
            % make the fitobject
            if ds.flags.fitRates
                ds.likeOpt(si).likeFun = ...
                    fittype('1 ./ sqrt(rmax .* (x .^ q ./ (x .^ q + c50 .^ q)) + rmin)', ...
                    'Coefficients', {'q', 'c50', 'rmin', 'rmax'});
            else
                ds.likeOpt(si).likeFun = ...
                    fittype('1 ./ sqrt(rmax .* (x .^ q ./ (x .^ q + c50 .^ q)) + rmin)', ...
                    'Coefficients', {'q', 'c50'}, 'Problem', {'rmin', 'rmax'});
            end
        else
            % in this case, use fmincon to attempt to minimize
            % least-squared error between the fitted values and functional
            % form
            
            % options to fmincon
            ds.likeOpt(si).options = ...
                optimset('Algorithm', 'active-set', 'Display', 'off', ...
                'MaxFunEvals', 5000);
            
            % define functional form
            if ds.flags.fitRates
                ds.likeOpt(si).likeFun = ...
                    @(c, q, c50, rMin, rMax) ...
                    (1 ./ sqrt(rMax .* (c .^ q ./ (c .^ q + c50 .^ q + rMin))));
            else
                ds.likeOpt(si).likeFun = ...
                    @(c, q, c50) ...
                    (1 ./ sqrt(ds.likeOpt(si).R(2) .* (c .^ q ./ ...
                    (c .^ q + c50 .^ q)) + ds.likeOpt(si).R(1)));
            end

			%%% I have not actually included this yet.
        end
    end
end

%% organize
ds = orderDataStruct(ds);
fprintf('done.\n');
