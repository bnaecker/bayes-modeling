function ds = plotPsychometricFunctions(ds)
%
% FUNCTION ds = plotPsychometricFunctions(ds)
%
% Plots the psychometric functions for the Bayesian ideal observer model,
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
    
    % preallocate cell arrays to hold the axes handles
    ds.handles(si).pmfAx = ds.pmfs(si).proportionTestFaster;
    
    % loop over each unique triplet of conditions, starting with reference
    % contrast (a loop over figures)
    for refCNum = 1:ds.info.nUniqueRefContrasts
        % start the index into the subplots
        spi = 0;
        
        % create a figure for each subject and reference contrast
        ds.handles(si).pmfFig(refCNum) = figure;
        set(ds.handles(si).pmfFig(refCNum), 'Color', 'w', 'NumberTitle', ...
            'off', 'Name', ['subject ' num2str(si) ' PMFs, reference ' ...
            'contrast = ' num2str(ds.data(si).refContrasts(refCNum))], ...
            'Units', 'normalized');
        
        % create labels to indicate direction of increasing speed and test
        % contrast
        ds.handles(si).pmfTextAx(refCNum) = axes('Position', ...
            [0.05 0.95 0.05 0.05], 'Units', 'Normalized', 'Visible', 'off');
        hold on;
        ds.handles(si).pmfRefVelText(refCNum) = text(0, 0.2, 'ref. vel.', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).pmfRefVelArrow(refCNum) = text(0.5, 0, '\downarrow', ...
            'FontSize', 12, 'FontWeight', 'bold');
        ds.handles(si).pmfTestCArrow(refCNum) = text(0, 0.6, ['test con. ' ...
            '\rightarrow'], 'FontSize', 12, 'FontWeight', 'bold');
        
        % loop over reference velocity
        for refVNum = 1:ds.info.nUniqueRefVels
            % loop over test contrast
            for testCNum = 1:ds.info.nUniqueContrasts
                
                % indices into this triplet condition
                refCInds = ds.data(si).refC == ...
                    ds.data(si).refContrasts(refCNum);
                refVInds = ds.data(si).refV == ...
                    ds.data(si).refVels(refVNum);
                testCInds = ds.data(si).testC == ...
                    ds.data(si).testContrasts(testCNum);
                
                % get the test and reference velocities for this triplet
                % condition, and transform them from the logarithmic domain
                testVels = unique(ds.data(si).testV(...
                        refCInds & refVInds & testCInds));
                refVels = ds.data(si).refVels; 
                if ds.info.isSimData
                    transTestVels = testVels;
                    transRefVels = refVels;
                else
                    transTestVels = ds.velTrans.iTransFun(testVels);
                    transRefVels = ds.velTrans.iTransFun(refVels);
                end
                
                % check that there is data for this condition
                if ~isempty(transTestVels)
                    
                    % preallocate the array to hold the proportion of times
                    % the test grating was reported as moving faster
                    ds.pmfs(si).proportionTestFaster{refCNum, refVNum, testCNum} = ...
                        nan(1, length(transTestVels));
                    
                    % grab the speed and choice data explicitly, so that we 
                    % can do maximum-likelihood fitting of the pmf
                    if ds.info.isSimData
                        speed = ds.data(si).testV(refCInds & refVInds & testCInds);
                    else
                        speed = ds.velTrans.iTransFun(...
                            ds.data(si).testV(refCInds & refVInds & testCInds));
                    end
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
                        startPoint = [transRefVels(refVNum); ...
                            1 / ds.data(si).testContrasts(testCNum)];
                        
                        % fit options
                        fitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                            'StartPoint', startPoint, 'Lower', [0 0]);
                        
                        % create the fit object
                        fitObj = fittype(weibullFun, 'Independent', 'x', ...
                            'Coefficients', {'lambda', 'k'}, 'Options', fitOpt);
                        
                        % use the curve fitting toolbox to create the fit
                        try
                        theFit = fit(transTestVels, ...
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
                            linspace(0.5 * transRefVels(refVNum) , ...
                            2 * transRefVels(refVNum), 1000);
                        
                        % use the fit object to create the pmf
                        ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum} = ...
                            theFit(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum});
                    else
                        % no curve fitting toolbox, so we'll fit the pmf by
                        % maximum-likelihood
                        
                        % initial parameter values, with bounds
                        prs0 = [transRefVels(refVNum); ...
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
                            linspace(0.5 * transRefVels(refVNum), ...
                            2 * transRefVels(refVNum), 1000);
                        
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
                        x1 = transRefVels(refVNum);
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
                        x1 = transRefVels(refVNum);
                        a1 = ds.params(si).slopeHat(refVNum, 1);
                        sig1 = ds.params(si).likeWidth(...
                            ds.data(si).refContrasts(refCNum) == ...
                            ds.data(si).testContrasts, 1);
                        sig2 = ds.params(si).likeWidth(testCNum, 1);
                        
                        % then for the test stimulus
                        x2 = ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum};
                        a2 = interp1(transRefVels, ...
                            ds.params(si).slopeHat(:, 1), x1, 'linear', 'extrap');
                        
                        % compute the predicted psychometric function
                        ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum} = ...
                            sort(normcdf(((x2 - x1 + a1 .* sig1 .^ 2 - ...
                            a2 .* sig2 .^ 2) ./ sqrt(sig1 .^ 2 + sig2 .^ 2))));
                    else
                        % optional functionality
                    end
                end
                
                % plot the psychometric functions, if requested
                if ds.flags.plotPMFs
                    
                    % index into the correct subplot
                    spi = spi + 1;
                    ds.handles(si).pmfAx{refCNum, refVNum, testCNum} = ...
                        subplot(ds.info.nUniqueRefVels, ds.info.nUniqueContrasts, ...
                        spi);
                    hold on;
                    
                    % label the rows and columns correctly (only on the
                    % outside)
                    if mod(spi, ds.info.nUniqueContrasts) == 1
                        % the subplot is along the left column
                        ds.handles(si).pmfVLabels(refCNum, refVNum) = ...
                            ylabel({[num2str(transRefVels(refVNum)) '       '], ...
                            ' (deg s^{-1})          '}, ...
                            'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 0);
                    elseif any(spi == (ds.info.nUniqueRefVels - 1) * ...
                            ds.info.nUniqueContrasts + 1 : ...
                            ds.info.nUniqueRefVels * ds.info.nUniqueContrasts)
                        
                        % the subplot is along the bottom row
                        ds.handles(si).pmfCLabels(refCNum, refVNum) = ...
                            xlabel(num2str(ds.data(si).testContrasts(testCNum)), ...
                            'FontSize', 12, 'FontWeight', 'bold');
                    end
                    
                    % fill the subplots
                    if ~isempty(testVels)
                        
                        % set the axis limits
                        set(ds.handles(si).pmfAx{refCNum, refVNum, testCNum}, ...
                            'TickDir', 'out', 'XLim', ...
                            [ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}(1) ...
                            ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}(end)], ...
                            'Units', 'normalized'); axis square; box on;
                        
                        
                        % plot the fitted pmf
                        ds.handles(si).pmfDataFit(refCNum, refVNum, testCNum) = ...
                            plot(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
                            ds.pmfs(si).dataPmfs{refCNum, refVNum, testCNum}, ...
                            'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none', ...
                            'Color', 'r');
                        
                        % plot the model-predicted pmf
                        ds.handles(si).pmfModelFit(refCNum, refVNum, testCNum) = ...
                            plot(ds.pmfs(si).pmfDomain{refCNum, refVNum, testCNum}, ...
                            ds.pmfs(si).modelPmfs{refCNum, refVNum, testCNum}, ...
                            'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none', ...
                            'Color', 'b');
                        
                        % plot the data
                        ds.handles(si).pmfDataLines(refCNum, refVNum, testCNum) = ...
                            plot(transTestVels', ...
                            ds.pmfs(si).proportionTestFaster{...
                            refCNum, refVNum, testCNum}(1, :), ...
                            'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', ...
                            'k', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
                        
                    else
                        % fill the axes with 'No Data' text
                        set(ds.handles(si).pmfAx{refCNum, refVNum, testCNum}, ...
                            'Visible', 'off');
                        text(0.5, 0.5, 'no data');
                    end
                    
                    % make one legend for each figure
                    if spi == ds.info.nUniqueContrasts
                        legend({'wbl fit', 'model pmf', 'data'}, ...
                            'Location', 'NorthEastOutside');
                    end
                end
            end
        end
    end
end