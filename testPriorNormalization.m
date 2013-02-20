%%
xPri = var.data(1).refVels;
ySlopes = var.params(1).slopeHat(1,:);

v0 = .3;
fnlin = @(v)(log(1+v/v0));
finv = @(u)(v0*(exp(u)-1));

vv = 0:.1:12;
plot(vv, finv(fnlin(vv)), vv, vv, 'k-o');


%%  Make prior sum to 1 in linear velocity

xvel = finv(xPri);  % x axis back in linear units

dx = .01;
xx = min(xvel):dx:max(xvel);
%xx = 2*dx:dx:max(xvel);
yysl = interp1(xvel,ySlopes,xx,'cubic', 'extrap');

ycum = cumsum(yysl)*dx;
ypri = exp(-ycum);
ypri = ypri./sum(ypri)/dx;

plot(xx, ypri);
loglog(xx, ypri);


%%  Make prior sum to 1 in log velocity

xvel = xPri;  % x axis still in log units

dx = .01;
xx = min(xvel):dx:max(xvel);
xx = finv(.4):dx:(max(xvel)+dx);
yysl = interp1(xvel,ySlopes,xx,'spline', 'extrap');

ycum = cumsum(yysl)*dx;
ypri = exp(-ycum);
ypri = ypri./sum(ypri)/dx;

plot(xx, ypri);
ypriVals = interp1(xx,ypri,xvel,'linear');

ii = (xx>min(xvel)-.2);
loglog(finv(xx(ii)), ypri(ii), finv(xvel), ypriVals, 'bo');



