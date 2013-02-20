function var = setupOptimization(var)
%
% Setup log-likelihood optimization
%
% (c) benjamin naecker UT Austin 20 May 2011 benjamin.naecker@gmail.com

for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    %% Setup optimization for negative log-likelihood of data
    % Log-likelihood function handle
    if var.flags.fitG
    else
        var.llOpt(subject).llFun = @(prs) (-logli_var(prs, var));
    end
    
    % Options to fmincon
    if var.info.nBoots == 1
        var.llOpt(subject).options = optimset('Algorithm', 'active-set', 'Display', 'iter', ...
            'MaxFunEvals', 5000, 'MaxIter', 500);
    else
        var.llOpt(subject).options = optimset('Algorithm', 'active-set', 'Display', 'off', ...
            'MaxFunEvals', 20000, 'MaxIter', 1000);
    end
    
    % Number of allowed failures-to-converge before the fitting asks you if
    % you'd like to continue
    var.llOpt(subject).numAllowedFails = 1;
    
    % Guesses at likelihood widths and priors
    hGuess = 1 .* ones(var.info.nUniqueContrasts, 1);
    hLB = .01 .* ones(var.info.nUniqueContrasts, 1);
    hUB = 100 .* ones(var.info.nUniqueContrasts, 1);
    priorGuess = .5 .* ones(var.info.nUniqueRefVels, 1);
    priorLB = -100 .* ones(var.info.nUniqueRefVels, 1);
    priorUB = 100 .* ones(var.info.nUniqueRefVels, 1);
    if var.flags.fitG
        % Add guesses and bounds for g(v)
    end
    
    % Organize into prs, LB, UB
    var.llOpt(subject).prs0 = [hGuess; priorGuess]; 
    var.llOpt(subject).LB = [hLB; priorLB];
    var.llOpt(subject).UB = [hUB; priorUB];
    if var.flags.fitG
        % Append guesses and bounds for g(v)
    end
    
    %% Setup optimization for LSE fit to h(c)
    if var.flags.fitLikeFun
        % Options to fmincon
        var.hOpt(subject).options = optimset('Algorithm', 'active-set', 'Display', 'off', 'MaxFunEvals', 500);
        
        if var.flags.fitRates
            % Guesses at parameters [q; c50; rMin; rMax]
            var.hOpt(subject).prs0 =  [2; .2; 10; 50;];
            var.hOpt(subject).LB = [1e-15; .01; 1; 20];
            var.hOpt(subject).UB = [100; 1; 20; 200];
            
            % Functional form
            var.hOpt(subject).hFun = @(c, q, c50, rMin, rMax) (1 ./ sqrt( rMax .* (c .^ q ./ (c .^ q + c50 .^ q + rMin))) );
        else
            % Guesses at parameters [q; c50]
            var.hOpt(subject).prs0 = [2; .2];
            var.hOpt(subject).LB = [1e-15; 1e-15];
            var.hOpt(subject).UB = [100; 1];
            
            % Functional form
            var.hOpt(subject).R = [.1; 30];
            var.hOpt(subject).hFun = @(c, q, c50) (1 ./ sqrt( (var.hOpt(subject).R(2) .* (c .^ q ./ (c .^ q + c50 .^ q))) + var.hOpt(subject).R(1)) );
        end
    end
end

%% Organize
var = organizeVar(var);