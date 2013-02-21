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
    elseif strcmpi(ds.flags.priorType, 'loglinear')
        ds.llOpt(si).llFun = @(prs) (-loglike_loglinear(prs, ds));
    elseif strcmpi(ds.flags.priorType, 'gaussian')
        ds.llOpt(si).llFun = @(prs) (-loglike_gaussian(prs, ds));
    else
        error('runSpeedDiscriminationModel:setupOptimization:unknownPrior', ...
            'The supported prior types are "loglinear" or "gaussian"');
    end
    
    %% options to fmincon
    if ds.info.nBoots == 1
        ds.llOpt(si).options = ...
            optimset('Algorithm', 'active-set', 'Display', 'Iter', ...
            'MaxFunEvals', 5000, 'MaxIter', 500);
    else
        ds.llOpt(si).options = ...
            optimset('Algorithm', 'active-set', 'Display', 'off', ...
            'MaxFunEvals', 5000, 'MaxIter', 500);
    end
    
    %% initial values for likelihood widths and prior parameters
    if strcmpi(ds.flags.priorType, 'loglinear')
        % guesses and bounds to likelihood widths
        % [1.3850 1.3290 1.2770 1.1118 0.9209 0.8631 0.7573]'		%% debugging an old simulation
        likeGuess =  1 .* ones(ds.info.nUniqueContrasts, 1);
        likeLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        likeUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        
        % guesses and bounds to prior slopes
        % [6.5 8.5 9 7.9 2 0.05]'		%% debugging an old simulation
        priorGuess = 2 .* ones(ds.info.nUniqueRefVels, 1);
        priorLB = -100 .* ones(ds.info.nUniqueRefVels, 1);
        priorUB = 100 .* ones(ds.info.nUniqueRefVels, 1);
        
    elseif strcmpi(ds.flags.priorType, 'gaussian')
        % guesses and bounds to likelihood widths
        likeGuess = 1 .* ones(ds.info.nUniqueContrasts, 1);
        likeLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        likeUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        
        % guesses and bounds to prior variance
        priorGuess = 10;
        priorLB = 1e-5;
        priorUB = 1e5;
        
    else
        % optional functionality
    end
    
    % organize the guesses
    ds.llOpt(si).prs0 = [likeGuess; priorGuess];
    ds.llOpt(si).LB = [likeLB; priorLB];
    ds.llOpt(si).UB = [likeUB; priorUB];
    
    
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
