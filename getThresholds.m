function var = getThresholds(var)
%
% var = getThresholds(var) calculates absolute and relative thresholds
% (Weber fractions) for subjects performing the speed discrimination task
% in Stocker & Simoncelli 2006
%
% (c) benjamin naecker UT Austin 31 May 2011 benjamin.naecker@gmail.com

for subject = 1:var.info.nSubjects
    var.info.currentSubject = subject;

 %% Plot PMFs if requested
 % Setup cell array to receive variable proportions
 var.stats(subject).testSeenFaster = cell(var.info.nUniqueRefConts,...
     var.info.nUniqueContrasts, var.info.nUniqueRefVels);
 
% Loop through each reference contrast
 for refContrastNum = 1:var.info.nUniqueRefConts
     % Create figure for each subject and reference contrast
     o.handles(subject).pmfFig= figure;
     set(gcf, 'Name', ['Subject ' num2str(subject) ': PMFs, Ref Contrast = '...
         num2str(var.data(subject).refContrasts(refContrastNum))], 'Color', 'w',...
         'Units', 'Normalized');
     
     % Create labels for direction of increasing ref speed and test
     % contrast
     o.handles(subject).pmfTextAx = axes('Position', [.05 .05 .05 .05], ...
         'Units', 'Normalized', 'Visible', 'off'); hold on;
     o.handles(subject).pmfRefVelArrow = text(.5, .5, 'Ref Vel. \rightarrow',...
         'FontSize', 12, 'FontName', 'Helvetica');
     o.handles(subject).pmfTestCArrow = text(0, 1, 'Test Contr. \rightarrow',...
         'FontSize', 12, 'FontName', 'Helvetica', 'Rotation', 90);
     
     % Set index to subplots
     subIndex = 0;
     
     % Go through each condition (ref vel & test contr), and plot the
     % number of times the subject indicated the "reference" grating was
     % seen tomove faster as a function of the test grating speed (i.e., a
     % psychometric function)
     for testContrastNum = 1:var.info.nUniqueContrasts
         for refVelocityNum = 1:var.info.nUniqueRefVels
             % Get test velocities at this condition
             testVels = var.data(subject).testV(...
                            var.data(subject).refC == ...
                                var.data(subject).refContrasts(refContrastNum) & ...
                            var.data(subject).testC == ...
                                var.data(subject).testContrasts(testContrastNum) & ...
                            var.data(subject).refV == ...
                                var.data(subject).refVels(refVelocityNum) ,:);
             
             % At each test velocity, calculate p(test seen faster)
             var.stats(subject).testSeenFaster{refContrastNum, testContrastNum, ...
                 refVelocityNum}(testVelocityNum)
         end
     end
 end
 
end