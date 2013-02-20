function var = organizeVar(var)
%
% var = organizeVar(var) is a helper function that orders all the structs
% and substructs of the variable 'var', containing all the information for
% the analysis of speed discrimination data
%
% (c) benjamin naecker UT Austin 2 Jun 2011 benjamin.naecker@gmail.com

var = orderfields(var);
names = fieldnames(var);
for fieldNum = 1:length(names)
    if isstruct(var.(names{fieldNum}))
        var.(names{fieldNum}) = orderfields(var.(names{fieldNum}));
    end
end