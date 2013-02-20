function ds = orderDataStruct(ds)
%
% Helper function to keep the data structure organized
%
% (c) bnaecker@stanford.edu 17 Sep 2012

ds = orderfields(ds);
names = fieldnames(ds);
for fieldNum = 1:length(names)
    if isstruct(ds.(names{fieldNum}))
        ds.(names{fieldNum}) = orderfields(ds.(names{fieldNum}));
    end
end