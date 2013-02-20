function ds = fitLikelihoodFunction(ds)
%
% FUNCTION ds = fitLikelihoodFunction(ds)
%
% Fits the requested functional form to the likelihood widths by minimizing
% the least-squared error. The function itself is taken from Sclar,
% Maunsell, and Lennie, 1990. 
%
% (c) bnaecker@stanford.edu 11 Nov 2012

%% check that the user requested a functional fit
if ~ds.flags.fitLikeFun
    return
end

%% notify
fprintf('Fitting functional form to contrast-dependence of likelihood widths ... ');

%% fit
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    % try to use curve fitting toolbox
    if ds.flags.useCFit
        % use fit.m
        if ds.flags.fitRates
            ds.likeOpt(si).fitObj = fit(ds.data(si).testContrasts, ...
                ds.params(si).likeWidth(:, 1), ...
                ds.likeOpt(si).likeFun, ds.likeOpt(si).fitOptions);
        else
            ds.likeOpt(si).fitObj = fit(ds.data(si).testContrasts, ...
                ds.params(si).likeWidth(:, 1), ...
                ds.likeOpt(si).likeFun, ds.likeOpt(si).fitOptions, ...
                'problem', {ds.likeOpt(si).R(2), ds.likeOpt(si).R(1)});
        end
        
        % define domain of function (used to plot)
        ds.likeOpt(si).domain = linspace(min(ds.data(si).testContrasts), ...
            max(ds.data(si).testContrasts), 1000);
    else
        
        % setup objective function
        ds.likeOpt(si).likeObjFun = @(prs) squaredError(prs, ds);
        
        % minimized squared error with fmincon
        ds.likeOpt(si).prsHat = fmincon(ds.likeOpt(si).likeObjFun, ...
            ds.likeOpt(si).prs0, [], [], [], [], ...
            ds.likeOpt(si).LB, ds.likeOpt(si).UB, [], ...
            ds.likeOpt(si).options);
        
        % save output
        ds.params(si).q = ds.likeOpt(si).prsHat(1);
        ds.params(si).c50 = ds.likeOpt(si).prsHat(2);
        if ds.flags.fitRates
            ds.params(si).rMin = ds.likeOpt(si).prsHat(3);
            ds.params(si).rMax = ds.likeOpt(si).prsHat(4);
        end
        
        % define domain of the function
        ds.likeOpt(si).domain = linspace(min(ds.data(si).testContrasts), ...
            max(ds.data(si).testContrasts), 1000);
        
        % values of the likelihood function at the actual contrasts
        if ds.flags.fitRates
            ds.data(si).likeFunVals = ds.likeOpt(si).likeFun(ds.likeOpt(si).domain, ...
                ds.params(si).q, ds.params(si).c50, ds.params(si).rMin, ...
                ds.params(si).rMax);
        else
            ds.data(si).likeFunVals = ds.likeOpt(si).likeFun(ds.likeOpt(si).domain, ...
                ds.params(si).q, ds.params(si).c50) ./ ds.data(si).G;
        end
    end
end

%% organize
ds = orderDataStruct(ds);
fprintf('done.\n');