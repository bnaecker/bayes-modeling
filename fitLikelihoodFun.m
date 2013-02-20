function ds = fitLikelihoodFun(ds)
%
% function ds = fitLikelihoodFun(ds)
%
% Fits a function defined in setupData.m to the likelihood width by
% minimizing least-squared error. The function is taken from Sclar,
% Maunsell, and Lennie, 1990.
%
% (c) bnaecker@stanford.edu 24 Mar 2012

%% Check
if ~ds.flags.fitLikeFun
    return
end

%% Notify
fprintf('\nFitting functional form to h(c)...');

%% Fit
for subject = 1:ds.info.nSubjects
    ds.info.currentSubject = subject;
    
    %% Try to fit with Curve Fitting Toolbox
    v = ver('curvefit');
    if ~isempty(v)
        if ds.flags.fitRates
            ds.hOpt(subject).fitObj = fit(ds.data(subject).testContrasts, ...
                ds.params(subject).likeWidth(1, :)', ...
                ds.hOpt(subject).hFun, ds.hOpt(subject).fitOptions);
        else
            ds.hOpt(subject).fitObj = fit(ds.data(subject).testContrasts, ...
                ds.params(subject).likeWidth(1, :)', ...
                ds.hOpt(subject).hFun, ds.hOpt(subject).fitOptions, ...
                'problem', {ds.hOpt(subject).R(2), ds.hOpt(subject).R(1)});
        end
        
        %% Define domain of function (used to plot)
        ds.hOpt(subject).domain = linspace(min(ds.data(subject).testContrasts), ...
            max(ds.data(subject).testContrasts), 1000);
        
    else % No cfit toolbox, use fmincon
        %% Setup objective function with most current version of data structure
        ds.hOpt(subject).hObjFun = @(prs) leastSquaredErr(prs, ds);
        
        %% Minimized LSE with fmincon
        ds.hOpt(subject).prsHat = fmincon(ds.hOpt(subject).hObjFun, ...
            ds.hOpt(subject).prs0, [], [], [], [], ...
            ds.hOpt(subject).LB, ds.hOpt(subject).UB, [], ...
            ds.hOpt(subject).options);
        
        %% Save output
        ds.params(subject).q = ds.hOpt(subject).prsHat(1);
        ds.params(subject).c50 = ds.hOpt(subject).prsHat(2);
        if ds.flags.fitRates
            ds.params(subject).rMin = ds.hOpt(subject).prsHat(3);
            ds.params(subject).rMax = ds.hOpt(subject).prsHat(4);
        end
        
        %% Define domain of function
        ds.hOpt(subject).domain = linspace(min(ds.data(subject).testContrasts), ...
            max(ds.data(subject).testContrasts), 1000);
        
        %% Values of h(c)
        if ds.flags.fitRates
            ds.data(subject).hFunVals = ds.hOpt(subject).hFun(ds.hOpt(subject).domain, ...
                ds.params(subject).q, ds.params(subject).c50, ds.params(subject).rMin, ...
                ds.params(subject).rMax);
        else
            ds.data(subject).hFunVals = ds.hOpt(subject).hFun(ds.hOpt(subject).domain, ...
                ds.params(subject).q, ds.params(subject).c50) ./ ds.data(subject).G;
        end
    end
end

%% Organize
ds = orderDataStruct(ds);

%% Notify
fprintf('done.\n');