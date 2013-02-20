function ds = setupOptimization(ds)
%
% FUNCTION ds = setupOptimization(ds)
%
% The function setupOptimization.m describes the optimization of the
% log-likelihood function, including starting values, bounds, and the form
% of the functions to be fitted.
%
% (c) bnaecker@stanford.edu 11 Feb 2012

fprintf('\nSetting up optimization parameters...');
for subject = 1:ds.info.nSubjects
    ds.info.currentSubject = subject;
    
    %% Setup maximization of the model itself
    % Log-likelihood function handle
    if ds.flags.fitG
    elseif strcmp(ds.flags.priorType, 'loglinear')
        ds.llOpt(subject).llFun = @(prs) (-logli_loglinear(prs, ds));
    elseif strcmp(ds.flags.priorType, 'gaussian')
        ds.llOpt(subject).llFun = @(prs) (-logli_gaussian(prs, ds));
    else
        error('RunSpeedDisriminationModel:setupOptimization:unknownPrior', ...
            'The supported prior types are "loglinear" or "gaussian"');
    end
    
    %% Options to fmincon
    if ds.info.nBoots == 1
        ds.llOpt(subject).options = ...
            optimset('Algorithm', 'active-set', 'Display', 'Iter', ...
            'MaxFunEvals', 5000, 'MaxIter', 500);
    else
        ds.llOpt(subject).options = ...
            optimset('Algorithm', 'active-set', 'Display', 'iter', ...
            'MaxFunEvals', 20000, 'MaxIter', 1000);
    end
    
    %% Guesses at likelihood widths and the prior parameters
    if strcmp(ds.flags.priorType, 'loglinear')
        % For the locally log-linear model, the prior is described by a
        % local slope around each reference velocity
        hGuess = 1 .* ones(ds.info.nUniqueContrasts, 1);
        hLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        hUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        priorGuess = .5 .* ones(ds.info.nUniqueRefVels, 1);
        priorLB = -100 .* ones(ds.info.nUniqueRefVels, 1);
        priorUB = 100 .* ones(ds.info.nUniqueRefVels, 1);
    else
        % For the full gaussian model, the prior is described by its
        % variance. The description of the likelihood function remains the
        % same.
        hGuess = 1 .* ones(ds.info.nUniqueContrasts, 1);
        hLB = .001 .* ones(ds.info.nUniqueContrasts, 1);
        hUB = 100 .* ones(ds.info.nUniqueContrasts, 1);
        priorGuess = 10;
        priorLB = 1e-5;
        priorUB = 1e5;
    end
    
    % Organize the guesses
    ds.llOpt(subject).prs0 = [hGuess; priorGuess];
    % DEBUGGING USING SIMULATION VALUEs
    %%%%% THESE WORK %%%%%
    % LOGLINEAR PRIOR
%     data.llOpt(subject).prs0 = [2; 1; .5; .25; .2; logspace(-.5, -1, 6)'];
    % GAUSSIAN PRIOR
%     data.llOpt(subject).prs0 = [3; 6; 5];
%     data.llOpt(subject).prs0 = [10; 5; 2.5; 1.25; 1; 5];
    ds.llOpt(subject).LB = [hLB; priorLB];
    ds.llOpt(subject).UB = [hUB; priorUB];
    
    %% Setup optimization of least-squares fit to h(c)
    if ds.flags.fitLikeFun
        % Everything needs these values
        if ds.flags.fitRates
            % Guesses at parameters [slope; c50; rMin; rMax]
            ds.hOpt(subject).prs0 = [2; .2; 10; 50];
            ds.hOpt(subject).LB = [1e-15; .01; 0; 1];
            ds.hOpt(subject).UB = [100; 1; 20; 200];
        else
            % Guesses at parameters [slope; c50]
            ds.hOpt(subject).prs0 = [2; .2];
            ds.hOpt(subject).LB = [1e-15; 1e-15];
            ds.hOpt(subject).UB = [100; 1];
            ds.hOpt(subject).R = [.1; 100];
        end
        
        % Try to use curve fitting toolbox
        v = ver('curvefit');
        if ~isempty(v)
            % fit options
            ds.hOpt(subject).fitOptions = ...
                fitoptions('Method', 'NonlinearLeastSquares',...
                'Lower', ds.hOpt(subject).LB, 'Upper', ds.hOpt(subject).UB, ...
                'StartPoint', ds.hOpt(subject).prs0);
            if ds.flags.fitRates
                % fit object itself
                ds.hOpt(subject).hFun = ...
                    fittype('1 ./ sqrt(rmax .* (x .^ q ./ (x .^ q + c50 .^ q)) + rmin)', ...
                    'Coefficients', {'q', 'c50', 'rmin', 'rmax'});
            else
                % fit object itself
                ds.hOpt(subject).hFun = ...
                    fittype('1 ./ sqrt(rmax .* (x .^ q ./ (x .^ q + c50 .^ q)) + rmin)', ...
                    'Coefficients', {'q', 'c50'}, 'Problem', {'rmax', 'rmin'});
            end
        else % we need to use fmincon to fit the likelihood function
            % Options to fmincon
            ds.hOpt(subject).options = ...
                optimset('Algorithm', 'active-set', 'Display', 'off', ...
                'MaxFunEvals', 500);
            
            if ds.flags.fitRates
                % Functional form
                ds.hOpt(subject).hFun = ...
                    @(c, q, c50, rMin, rMax) ...
                    ( 1 ./ sqrt( rMax .* ( c .^ 1 ./ (c .^ q + c50 .^ q + rMin))) );
            else
                % Functional form
                ds.hOpt(subject).hFun = ...
                    @(c, q, c50) ( 1 ./ sqrt( (ds.hOpt(subject).R(2) .* (c .^ q ./ ...
                    (c .^ q + c50 .^ q))) + ds.hOpt(subject).R(1)) );
            end
        end
    end
end

%% Organize
ds = orderDataStruct(ds);
fprintf('done.');