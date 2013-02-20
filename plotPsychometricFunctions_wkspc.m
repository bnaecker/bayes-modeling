function [data, x, y] = plotPsychometricFunctions(data, varargin)
%
% function data = plotPsychometricFunctions(data, varargin)
%
% Plot psychometric functions for the speed discrimination model.
%
% INPUT:    data - the standard data structure
%
%           varargin - string describing whether to show the Weibull fits
%           to data, the PMFs predicted from the model, or both. {'fits',
%           'model', 'both'}
%
% OUTPUT:   data - the standard data structure
%
% See also: RunSpeedDiscriminationModel.m, plotDistFits.m
%
% (c) bnaecker@stanford.edu 2 May 2012

%% Parse input
if isempty(varargin)
    warning('RunSpeedDiscriminationModel:plotPsychometricFunctions:tooManyInputArgs', ...
        'No extra arguments given, showing both Weibull fits and model-predicted PMFS');
    plotType = 'both';
elseif length(varargin) > 1
    error('RunSpeedDiscriminationModel:plotPsychometricFunctions:tooManyInputArgs', ...
        ['Extra arguments must be a single string defining whether to plot the ' ...
        'Weibull "fits", the "model" predicted PMFs, or "both"']);
else
    plotType = varargin{1};
end

%% Loop over subjects and calculate the desired PMFs
for si = 1:data.info.nSubjects
    data.info.currentSubject = si;
    
    % Preallocate cell arrays to hold the PMF for each 3-condition
    % Recall: a 3-condition is defined by (refVel, refC, testC)
    data.pmfs(si).propTestFaster = cell(data.info.nUniqueRefContrasts, ...
        data.info.nUniqueRefVels, data.info.nUniqueContrasts);
    data.pmfs(si).meanProps = data.pmfs(si).propTestFaster;
    data.pmfs(si).propSDs = data.pmfs(si).meanProps;
    data.pmfs(si).dataPmfs = data.pmfs(si).meanProps;
    data.pmfs(si).meanDataPmf = data.pmfs(si).meanProps;
    data.pmfs(si).modelPmfs = data.pmfs(si).meanProps;
    data.pmfs(si).meanModelPmf = data.pmfs(si).meanProps;
    data.pmfs(si).pmfDomain = data.pmfs(si).meanProps;
    
    % Loop through each reference contrast
    for refContrastNum = 1:data.info.nUniqueRefContrasts
        if data.flags.plotPMFs
            % Create a figure for each subject and reference contrast
            data.handles(si).pmfFig = figure;
            set(gcf, 'Name', ['Subject = ' num2str(si) ' PMFs, Ref Contrasts = ' ...
                num2str(data.data(si).refContrasts(refContrastNum))], 'Color', 'w', ...
                'Units', 'Normalized', 'NumberTitle', 'off');
            
            % Create labels for direction of increasing ref speed and test
            % contras
            data.handles(si).pmfTextAx = axes('Position', [.05 .95 .05 .05], ...
                'Units', 'Normalized', 'Visible', 'off'); hold on;
            data.handles(si).pmfRefVelTxt = text(0, .2, 'Ref Vel.', ...
                'FontSize', 12, 'FontName', 'Helvetica');
            data.handles(si).pmfRefVelArrow = text(.5, 0, '\downarrow', ...
                'FontSize', 12, 'FontName', 'Helvetica');
            data.handles(si).pmfTestCArrow = text(0, .6, 'Test Contr. \rightarrow', ...
                'FontSize', 12, 'FontName', 'Helvetica');
        end
        
        % Set index to subplots
        subIndex = 0;
        
        % Go through each condition (ref vel & test contr), and plot the
        % number of times the subject indicated the "test" grating was
        % seen to move faster as a function of the test grating speed (i.e., a
        % psychometric function)
        for refVelocityNum = 1:data.info.nUniqueRefVels
            for testContrastNum = 1:data.info.nUniqueContrasts
                
                % Indices to this 3-condition
                refCInds = data.data(si).refC == data.data(si).refContrasts(refContrastNum);
                refVInds = data.data(si).refV == data.data(si).refVels(refVelocityNum);
                testCInds = data.data(si).testC == data.data(si).testContrasts(testContrastNum);
                
                % Get test velocities at this 3-condition
                testVels = unique(data.data(si).testV(refCInds & refVInds & testCInds, :));
                
                % Check that we have any data
                if ~isempty(testVels)
                    
                    %% The fits to data
                    if any(strcmp(plotType, {'fits', 'both'}))
                        % Preallocate the bootstrap array
                        data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, ...
                            testContrastNum} = nan(data.info.nBoots, length(testVels));
                        data.pmfs(si).pmfPrs{refContrastNum, refVelocityNum, ...
                            testContrastNum} = nan(data.info.nBoots, 2);
                        data.pmfs(si).dataPmfs{refContrastNum, refVelocityNum, ...
                            testContrastNum} = nan(data.info.nBoots, 1000);
                        data.pmfs(si).modelPmfs{refContrastNum, refVelocityNum, ...
                            testContrastNum} = nan(data.info.nBoots, 1000);
                        
                        % Bootstrap an entire range of the test velocities, so
                        % that we can fit a PMF to it at the same time
                        for boot = 1:data.info.nBoots
                            for testVelocityNum = 1:length(testVels)
                                % Resample indices
                                INDS = refCInds & refVInds & testCInds & data.data(si).testV == ...
                                    testVels(testVelocityNum);
                                
                                % Calculate P(test seen faster) for each
                                % testVel. Choice == 1 means reference was
                                % faster, so we need choice == 0
                                data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, ...
                                    testContrastNum}(boot, testVelocityNum) = ...
                                    sum(data.data(si).T(INDS, boot)) ./ ...
                                    length(data.data(si).T(INDS, boot));
                            end
                            
                            % Fit requested cumulative distribution
                            % function to each bootstrapped set of testVels
                            prs0 = [data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum)); ...
                                1 / data.data(si).testContrasts(testContrastNum)];
                            LB = [1e-15; 1e-15];
                            UB = [100; 100];
                            
                            % Define the objective function for fmincon.
                            % It's the squared-error bewteen the CDF and
                            % the data.
                            fitFun = @(prs) fitCDF(prs, data.flags.pmfFitType, ...
                                data.vTrans.iTransFun(testVels'), data.pmfs(si).propTestFaster{ ...
                                refContrastNum, refVelocityNum, testContrastNum}(boot, :));
                            data.pmfs.pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, :) = ...
                                fmincon(fitFun, prs0, [], [], [], [], LB, UB, [], data.llOpt(si).options);
                            data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum} = ...
                                linspace(.5 * data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum)), ...
                                2 * data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum)), 1000);
                            data.pmfs(si).dataPmfs{refContrastNum, refVelocityNum, testContrastNum}(boot, :) = ...
                                cdf(data.flags.pmfFitType, data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}, ...
                                data.pmfs(si).pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, 1), ...
                                data.pmfs(si).pmfPrs{refContrastNum, refVelocityNum, testContrastNum}(boot, 2));
                            
                            if refContrastNum == 1 && refVelocityNum == 3 && testContrastNum == 3
                                x = testVels;
                                y = data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}(1, :);
                            end
                        end
                        
                        % Calculate mean/SDs of all the bootstrapped choice
                        % sets
                        data.pmfs(si).meanProps{refContrastNum, refVelocityNum, testContrastNum} = ...
                            mean(data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}, 1);
                        data.pmfs(si).propSDs{refContrastNum, refVelocityNum, testContrastNum} = ...
                            std(data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}, 1);
                        
                        % Calculate mean/SDs of bootstrapped PMF fits
                        data.pmfs(si).meanDataPmf{refContrastNum, refVelocityNum, testContrastNum} = ...
                            mean(data.pmfs(si).dataPmfs{refContrastNum, refVelocityNum, testContrastNum}, 1);
                    end
                    
                    %% The model-predicted PMFs
                    if any(strcmp(plotType, {'model', 'both'}))
                        if strcmp(data.flags.priorType, 'gaussian')
                        elseif strcmp(data.flags.priorType, 'loglinear')
                            % Pull out the slope of the prior and the width
                            % of the likelihood for this 3-condition
                            x1 = data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum));
                            a1 = mean(data.params(si).slopeHat(:, refVelocityNum));
                            sig1 = mean(data.params(si).likeWidth(:, ...
                                data.data(si).refContrasts(refContrastNum) == data.data(si).testContrasts));
                            sig2 = mean(data.params(si).likeWidth(:, testContrastNum));
                            
                            % For the test stimulus, the speed, and
                            % therefore slope of the prior, must be
                            % interpolated across the desired range
                            if isempty(data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum})
                                data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum} = ...
                                    linspace(.5 * data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum)), ...
                                    2 * data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum)), 1000);
                            end
                            x2 = data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum};
                            a2 = interp1(data.vTrans.iTransFun(data.data(si).refVels), ...
                                mean(data.params(si).slopeHat), x2, 'linear', 'extrap');
                            
                            % Compute the psychometric function
                            data.pmfs.modelPmfs{refContrastNum, refVelocityNum, testContrastNum} = ...
                                normcdf( sort((x2 - x1 + a1 .* sig1 .^ 2 - a2 .* sig2 .^ 2) ./ ...
                                sqrt(sig1 .^ 2 + sig2 .^ 2)));
                        else
                        end
                    end
                end
                
                %% Show the requested PMFs
                if data.flags.plotPMFs
                    % Index into the correct subplot
                    subIndex = subIndex + 1;
                    data.handles(si).pmfAx(refContrastNum, refVelocityNum, testContrastNum) = ...
                        subplot(data.info.nUniqueRefVels, data.info.nUniqueContrasts, subIndex); hold on;
                    
                    % Label the rows and columns correctly
                    if mod(subIndex, data.info.nUniqueContrasts) == 1
                        data.handles(si).pmfVLabels(refContrastNum, refVelocityNum) = ...
                            ylabel([num2str(data.vTrans.iTransFun(data.data(si).refVels(refVelocityNum))) ...
                            ' [deg s^{-1}]      '], 'FontSize', 12, 'FontWeight', 'Bold', 'Rotation', 0);
                        if subIndex == (data.info.nUniqueRefVels - 1) * data.info.nUniqueContrasts + 1
                            data.hadles(si).pmfCLabels(refContrastNum, testContrastNum) = xlabel( ...
                                num2str(data.data(si).testContrasts(testContrastNum)), ...
                                'FontSize', 12, 'FontWeight', 'Bold');
                        end
                    elseif any(subIndex == (data.info.nUniqueRefVels - 1) * data.info.nUniqueContrasts + 1 : ...
                            data.info.nUniqueRefVels * data.info.nUniqueContrasts)
                        data.handles(si).pmfVLabels(refContrastNum, testContrastNum) = xlabel(...
                            num2str(data.data(si).testContrasts(testContrastNum)), ...
                            'FontSize', 12, 'FontWeight', 'Bold');
                    end
                    if ~isempty(testVels)
                        set(data.handles(si).pmfAx(refContrastNum, refVelocityNum, testContrastNum), ...
                            'TickDir', 'out', 'XLim', [data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}(1) ...
                            data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}(end)], ...
                            'Units', 'Normalized');
                    end
                    
                    % Plot the data, plus the requested fits
                    if ~isempty(testVels)
                        data.handles(si).pmfData(refContrastNum, refVelocityNum, testContrastNum) = ...
                            plot(data.vTrans.iTransFun(testVels'), ...
                            data.pmfs(si).propTestFaster{refContrastNum, refVelocityNum, testContrastNum}(1, :), ...
                            'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 5);
                        if any(strcmp(plotType, {'fits', 'both'}))
                            data.handles(si).pmfDataFit(refContrastNum, refVelocityNum, testContrastNum) = ...
                                plot(data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}, ...
                                data.pmfs(si).dataPmfs{refContrastNum, refVelocityNum, testContrastNum}(1, :), ...
                                'LineStyle', '-.', 'LineWidth', 2, 'Marker', 'none', 'Color', [1 .3 .3]);
                        elseif any(strcmp(plotType, {'model', 'both'}))
                            data.handles(si).pmfModelFit(refContrastNum, refVelocityNum, testContrastNum) = ...
                                plot(data.pmfs(si).pmfDomain{refContrastNum, refVelocityNum, testContrastNum}, ...
                                data.pmfs(si).modelPmfs{refContrastNum, refVelocityNum, testContrastNum}(1, :), ...
                                'LineStyle', '-.', 'LineWidth', 2, 'Marker', 'none', 'Color', [.3 .3 1]);
                        end
                    else
                        set(data.handles(si).pmfAx(refContrastNum, refVelocityNum, testContrastNum), 'Visible', 'off');
                        data.handles(si).emptyText{refContrastNum, refVelocityNum, testContrastNum} = ...
                            text(.5, .5, 'No data');
                    end
                end
            end
        end
    end
end