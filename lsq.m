function err = lsq(prs, var)
%
% lsq calculates the squared-error between the functional form of the
% likelihood with h(c) using prs from fmincon and the actual data
%
% (c) benjamin naecker UT Austin 26 May 2011 benjamin.naecker@gmail.com

if var.flags.fitRates
    fit = var.hOpt(var.info.currentSubject).hFun(var.data(var.info.currentSubject).testContrasts, ...
        prs(1), prs(2), prs(3), prs(4))';
else
    fit = var.hOpt(var.info.currentSubject).hFun(var.data(var.info.currentSubject).testContrasts, ...
        prs(1), prs(2))';
end
    
data = var.params(var.info.currentSubject).hHat(1, :);

err = sum( (fit - data) .^ 2 );