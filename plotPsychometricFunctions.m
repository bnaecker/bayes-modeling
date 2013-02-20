function [ds, x, y] = plotPsychometricFunctions(ds)
%
% function ds = plotPsychometricFunctions(ds)
%
% Plots the psychometric functions for the speed discrimination model, 
% both as Weibull fits to the subjects' choice data and those predicted
% from the fits to the model itself.
%
% INPUT:    data - the standard data structure
%
% OUTPUT:   data - standard data structure
%
% See also: RunSpeedDiscriminationModel.m, plotDistFits.m
%
% (c) bnaecker@stanford.edu 18 May 2012

useCFit = 1;
notify = 0;

%% Loop over subjects to calculate the desired quantities for each. These 
%% include the proportion of times the test grating was reported to move 
%% faster (the PMF) and the fits to the data or from the model.
for si = 1:ds.info.nSubjects
    ds.info.currentSubject = si;
    
    % Preallocate cell array to hold the proportions, the fitted
    % parameters, and the PMFs themselves
    ds.pmfs(si).propTestFaster = cell(ds.info.nUniqueRefContrasts, ...
        ds.info.nUniqueRefVels, ds.info.nUniqueContrasts);
    ds.pmfs(si).dataPmfs = ds.pmfs(si).propTestFaster;
    ds.pmfs(si).modelPmfs = ds.pmfs(si).propTestFaster;
    ds.pmfs(si).pmfDomain = ds.pmfs(si).propTestFaster;
    
    % Loop through each unique 3-condition
    for refCNum = 1:ds.info.nUniqueRefContrasts
        % Start index into subplots
        subPlotIndex = 0;
        
        % Create a figure for each subject and reference contrast
        ds.handles(si).pmfFig = figure;
        set(gcf, 'Name', ['Subject = ' num2str(si) ' PMFs, Ref Contrasts = ' ...
            num2str(ds.data(si).refContrasts(refCNum))], 'Color', 'w', ...
            'Units', 'Normalized', 'NumberTitle', 'off');
        
        % Create labels for direction of increasing ref speed and test
        % contras
        ds.handles(si).pmfTextAx = axes('Position', [.05 .95 .05 .05], ...
            'Units', 'Normalized', 'Visible', 'off'); hold on;
        ds.handles(si).pmfRefVelTxt = text(0, .2, 'Ref Vel.', ...
            'FontSize', 12, 'FontName', 'Helvetica');
        ds.handles(si).pmfRefVelArrow = text(.5, 0, '\downarrow', ...
            'FontSize', 12, 'FontName', 'Helvetica');
        ds.handles(si).pmfTestCArrow = text(0, .6, 'Test Contr. \rightarrow', ...
            'FontSize', 12, 'FontName', 'Helvetica');
        
        for refVNum = 1:ds.info.nUniqueRefVels
            for testCNum = 1:ds.info.nUniqueContrasts
                % Indices to this condition
                refCInds = ds.data(si).refC == ...
                    ds.data(si).refContrasts(refCNum);
                refVInds = ds.data(si).refV == ...
                    ds.data(si).refVels(refVNum);
                testCInds = ds.data(si).testC == ...
                    ds.data(si).testContrasts(testCNum);
                
                % Get the test and reference velocities, and transform them
                % as needed
                testVels = unique(ds.data(si).testV(refCInds & refVInds & testCInds, :));
                refVels = ds.data(si).refVels;
                if ds.info.isSimData
                    transTestVels = testVels;
                    transRefVels = refVels;
                else
                    transTestVels = ds.vTrans.iTransFun(testVels);
                    transRefVels = ds.vTrans.iTransFun(refVels);
                end
                
                % Check that there is data for this condition
                if ~isempty(testVels)
                    % Preallocate this bootstrap array within the cell
                    ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum} = ...
                        nan(1, length(testVels));
                    
                    %% Must grab velocity and choice data explicitly for ML fitting of Weibulls
                    X = ds.data(si).testV(refCInds & refVInds & testCInds);
                    choice = ds.data(si).T(refCInds & refVInds & testCInds);
                    
                    % Loop over the test velocities
                    for testVNum = 1:length(testVels)
                        % Indices to the choice array
                        INDS = refCInds & refVInds & testCInds & ...
                            ds.data(si).testV == testVels(testVNum);
                        
                        %% Calculate proportion of stimuli in which the
                        % test grating was reported as faster. Choice
                        % == 1 means reference was faster, so look for
                        % choice == 0.
                        ds.pmfs(si).propTestFaster{refCNum, refVNum, ...
                            testCNum}(1, testVNum) = ...
                            sum(ds.data(si).T(INDS, 1)) / ...
                            length(ds.data(si).T(INDS, 1));
                    end
                    
                    if refCNum == 2 && refVNum == 1 && testCNum == 4
                        x = testVels;
                        y = ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum};
                    end
                    
                    %% Fit a cumulative Weibull to the choice data
                    if ~isempty(ver('curvefit')) && useCFit
                        
                        % Setup to the fit
                        weibull = '1 - exp(-(x ./ lambda) .^ k)';
                        startPoint = [transRefVels(refVNum); ...
                            1 / ds.data(si).testContrasts(testCNum)];
                        s = fitoptions('Method', 'NonlinearLeastSquares', ...
                            'StartPoint', startPoint, ...
                            'Lower', [0 0]);
                        fitObj = fittype(weibull, 'Independent', 'x', ...
                            'Coefficients', {'lambda', 'k'}, ...
                            'Options', s);
                        
                        % Do the fit from fitting toolbox (requires cfit
                        % toolbox)
                        theFit = fit(transTestVels, ...
                            ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum}', ...
                            fitObj);
                        
                        % Save the parameters
                        ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum} = ...
                            [theFit.k theFit.lambda];
                        
                        
                        ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum} = ...
                            linspace(.5 * transRefVels(refVNum), ...
                            2 * transRefVels(refVNum), 1000);
                        ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum} = ...
                            theFit(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum});
                    else
                        % Setup mimization parameters
                        prs0 = [transRefVels(refVNum); ...
                            ds.data(si).testContrasts(testCNum)];
                        LB = [1e-15; 1e-15];
                        UB = [100; 100];
                        
                        
                        % Define objective function, maximum-likelihood
                        fitFun = @(prs) fitCDF(prs, ds.flags.pmfFitType, ...
                            X, choice);
                        
                        % Fit the function
%                         prs0 = prs0 + (.1 .* prs0) .* rand(size(prs0));
                        [ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum}, residualSS] = ...
                            fmincon(fitFun, prs0, [], [], [], [], LB, UB, [], ds.llOpt(si).options);
                        
                        totalSS = sum( (ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum} - ...
                            mean(ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum})) .^ 2);
                        dataRSquare = 1 - (residualSS / totalSS);
                        
                        % Define the domain and the PMF itself
                        ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum} = ...
                            linspace(.5 * transRefVels(refVNum), ...
                            2 * transRefVels(refVNum), 1000);
                        ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum} = ...
                            cdf(ds.flags.pmfFitType, ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
                            ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum}(1), ...
                            ds.pmfs(si).pmfPrs{refCNum, refVNum, testCNum}(2));
                    end
                    
                    %% Compute the model-predicted PMF
                    if strcmp(ds.flags.priorType, 'gaussian')
                        % Pull out the width of the prior (same for all
                        % velocities) and the width of the likelihood for
                        % this condition, first for the reference stimulus
                        x1 = transRefVels(refVNum);
                        sig1 = ds.params(si).likeWidth(1, ...
                            ds.data(si).refContrasts(refCNum) == ds.data(si).testContrasts);
                        sig2 = ds.params(si).likeWidth(1, testCNum);
                        
                        % Then for the test stimulus
                        x2 = ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum};
                        
                        % Compute a1 and a2 (it's backwards from sigmas)
                        a2 = ds.params(si).gamma(1) .^ 2 + sig1 .^ 2;
                        a1 = ds.params(si).gamma(1) .^ 2 + sig2 .^ 2;
                        
                        % Compute the psychometric function
                        ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum} = ...
                            sort(normcdf((a2 .* x2 - a1 .* x1) ./ ...
                            sqrt(a2 .^ 2 .* sig2 .^ 2 + a2 .^ 2 .* sig1 .^ 2)));
%                         
%                         figure, subplot(211); plot(x2, a1); title('slope values');
%                             subplot(212), plot(data.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
%                                 data.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}); title('pmf');
%                             h = 1;
                    elseif strcmp(ds.flags.priorType, 'loglinear');
                        % Pull out the slope of the log-prior and the
                        % likelihood width for this condition, first
                        % for the reference stimulus
                        
                        x1 = transRefVels(refVNum);
                        a1 = ds.params(si).slopeHat(1, refVNum);
                        sig1 = ds.params(si).likeWidth(1, ...
                            ds.data(si).refContrasts(refCNum) == ds.data(si).testContrasts);
                        sig2 = ds.params(si).likeWidth(1, testCNum);
                        
                        % Then for the test stimulus
                        x2 = ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum};
                        a2 = interp1(transRefVels, ...
                            ds.params(si).slopeHat(1, :), x2, 'linear', 'extrap');
                        ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum} = ...
                            sort(normcdf( ((x2 - x1 + a1 .* sig1 .^ 2 - a2 .* sig2 .^ 2) ./ ...
                            sqrt(sig1 .^2 + sig2 .^ 2))));
                        
%                         if (refCNum == 2 && refVNum == 2 && testCNum == 4) || ...
%                                 (refCNum == 2 && refVNum == 1 && testCNum == 4)
%                             figure, subplot(211); plot(x2, a2); title('slope values');
%                             subplot(212), plot(data.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
%                                 data.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}); title('pmf');
%                             h = 1;
%                         end
                    else
                    end
                    
                    % Model coefficient of determination
                    if notify
                        f = interp1(x2, ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}, ...
                            transTtestVels, 'linear', 'extrap')';
                        residualSS = sum( (ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum} - f) .^ 2 );
                        totalSS = sum( (ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum} - ...
                            mean(ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum})) .^ 2 );
                        modelRSquare = 1 - (residualSS / totalSS);
                        fprintf('\nCondition (%d, %d, %d), Weibull R-squared : %.3f, Model R-squared: %.3f\n', ...
                            refCNum, refVNum, testCNum, dataRSquare, modelRSquare);
                    end
                end
                    
                %% Plot the pmfs, if requested
                if ds.flags.plotPMFs
                    
                    % Index into the correct subplot
                    subPlotIndex = subPlotIndex + 1;
                    ds.handles(si).pmfAx(refCNum, refVNum, testCNum) = ...
                        subplot(ds.info.nUniqueRefVels, ds.info.nUniqueContrasts, subPlotIndex);
                    hold on;
                    
                    % Label rows and columns correctly
                    if mod(subPlotIndex, ds.info.nUniqueContrasts) == 1
                        ds.handles(si).pmfVLabels(refCNum, refVNum) = ...
                            ylabel([num2str(transRefVels(refVNum)) ...
                            ' [deg s^{-1}]          '], 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 0);
                        if subPlotIndex == (ds.info.nUniqueRefVels - 1) * ds.info.nUniqueContrasts + 1
                            ds.hadles(si).pmfCLabels(refCNum, testCNum) = xlabel( ...
                                num2str(ds.data(si).testContrasts(testCNum)), ...
                                'FontSize', 12, 'FontWeight', 'Bold');
                        end
                    elseif any(subPlotIndex == (ds.info.nUniqueRefVels - 1) * ds.info.nUniqueContrasts + 1 : ...
                            ds.info.nUniqueRefVels * ds.info.nUniqueContrasts)
                        ds.handles(si).pmfVLabels(refCNum, testCNum) = xlabel(...
                            num2str(ds.data(si).testContrasts(testCNum)), ...
                            'FontSize', 12, 'FontWeight', 'Bold');
                    end
                    if ~isempty(testVels)
                        set(ds.handles(si).pmfAx(refCNum, refVNum, testCNum), ...
                            'TickDir', 'out', 'XLim', [ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}(1) ...
                            ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}(end)], ...
                            'Units', 'Normalized');
                    end
                    
                    % Plot the data and the fits
                    if ~isempty(testVels)
                        ds.handles(si).pmfData(refCNum, refVNum, testCNum) = ...
                            plot(transTestVels', ...
                            ds.pmfs(si).propTestFaster{refCNum, refVNum, testCNum}(1, :), ...
                            'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 5);
                        ds.handles(si).pmfDataFit(refCNum, refVNum, testCNum) = ...
                            plot(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
                            ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum}, ...
                            'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none', 'Color', [1 .3 .3]);
                        ds.handles(si).pmfModelFit(refCNum, refVNum, testCNum) = ...
                            plot(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
                            ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}, ...
                            'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none', 'Color', [.3 .3 1]);
                    else
                        set(ds.handles(si).pmfAx(refCNum, refVNum, testCNum), 'Visible', 'off');
                        ds.handles(si).emptyText{refCNum, refVNum, testCNum} = ...
                            text(.5, .5, 'No data');
                    end
                    
                    % Make one legend for each figure
                    if subPlotIndex == ds.info.nUniqueContrasts
                        legend({'Data', 'Weibull fit', 'Model'}, 'Location', 'NorthEastOutside');
                    end
                end
            end
        end
    end
end