function var = fitLikeFun(var)
%
% fitLikeFun fits a function defined in setupVar to the likelihood width
% using LSE. The function is taken from Sclar, Maunsell, and Lennie 1990
%
% (c) benjamin naecker UT Austin 26 May 2011 benjamin.naecker@gmail.com

%% Notify
fprintf('Fitting functional form to h(c)...');

%% Fit
for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;
    
    % Setup objective function with most current information
    var.hOpt(subject).hObjFun = @(prs)(lsq(prs, var));
    
    % Minimze with fmincon
    var.hOpt(subject).prsHat = fmincon(var.hOpt(subject).hObjFun, var.hOpt(subject).prs0, ...
        [], [], [], [], var.hOpt(subject).LB, var.hOpt(subject).UB, [], var.hOpt(subject).options);
    
    % Save output
    var.params(subject).q = var.hOpt(subject).prsHat(1);
    var.params(subject).c50 = var.hOpt(subject).prsHat(2);
    if var.flags.fitRates
        var.params(subject).rMin = var.hOpt(subject).prsHat(3);
        var.params(subject).rMax = var.hOpt(subject).prsHat(4);
    end
    
    % Define domain
    var.hOpt(subject).domain = linspace(min(var.data(subject).testContrasts),...
        max(var.data(subject).testContrasts), 1000);
    
    % Values of h(c)
    if var.flags.fitRates
        var.data(subject).hFunVals = var.hOpt(subject).hFun(var.hOpt(subject).domain,...
            var.params(subject).q, var.params(subject).c50, var.params(subject).rMin, var.params(subject).rMax) ...
            ./ var.data(subject).G;
    else
        var.data(subject).hFunVals = var.hOpt(subject).hFun(var.hOpt(subject).domain,...
            var.params(subject).q, var.params(subject).c50) ./ var.data(subject).G;
    end
end

%% Notify
fprintf('converged.\n');