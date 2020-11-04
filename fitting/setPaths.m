function ds = setPaths(ds)
%
% FUNCTION ds = setPaths(ds)
%
% The function setPaths sets the paths for the speed discrimination experiment.
% They are only loaded for the current MATLAB session, not saved.
%
% (c) bnaecker@stanford.edu 2013 
% 21 Jun 2013 - wrote it

% get the current directory, move up a few, and save
oldpwd = pwd;
cd('../../');
ds.info.baseDir = pwd;

% get the current matlab path info
p = regexp(path, ':', 'split');

% the subdirectories we care about
subdirs = cellfun(@(c) fullfile(ds.info.baseDir, c), ...
	{'code/fitting/', 'code/plotting/', 'code/testScripts/', 'data/'}, ...
	'UniformOutput', false);

% add the subdirs if they're not already there
for si = 1:length(subdirs)
	if all(cellfun(@(c) isempty(regexp(c, subdirs{si})), p, 'UniformOutput', true))
		addpath(subdirs{si});
	end
end

% change dirs back
cd(oldpwd);
