function err = leastSquaredErr(prs, ds)
%
% function err = leastSquaredErr(prs, data) calculates the squared error
% between the function form of the likelihood width h(c) using the
% parameters chosen by fmincon and the actual data
%
% (c) bnaecker@stanford.edu 24 Mar 2012

if ds.flags.fitRates
    fit = ds.hOpt(ds.info.currentSubject).hFun(...
        ds.data(ds.info.currentSubject).testContrasts, ...
        prs(1), prs(2), prs(3), prs(4))';
else
    fit = ds.hOpt(ds.info.currentSubject).hFun(...
        ds.data(ds.info.currentSubject).testContrasts, ...
        prs(1), prs(2))';
end

d = ds.params(ds.info.currentSubject).hHat(1, :);

err = sum( (fit - d) .^ 2 );