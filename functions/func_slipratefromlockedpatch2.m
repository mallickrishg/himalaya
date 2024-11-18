function [xgrid,sliprate] = func_slipratefromlockedpatch2(m,nmesh)
% compute fault creep rate as a function of horizontal location from the
% fault
% model parameters for inversion:
% m(1) - downdip locking limit in [km]
% m(2) - transition width (from lock to creep at plate rate) [km]
% m(3) - Vplate in [mm/yr]
% m(4) - fault dip in degrees
% m(5) - plate thickness in [km]
% AUTHOR:
% Rishav Mallick, JPL, 2024

% nmesh = 1000;

xlock = m(1);
w = m(2);
Vplate = m(3); % plate velocity in [mm/yr]
dip = m(4);

sliprate = zeros(nmesh,1);
zetamesh = linspace(0,300,nmesh)'; % this is distance along the fault
% project zeta -> x coordinate and convert to meters
xgrid = zetamesh.*cosd(dip).*1e3;

creepindex = (zetamesh > xlock) & zetamesh <= (xlock + w);
% using the traction free boundaries for the stress-driven creeping
% segments while the forced slip segments are either no slip or full
% plate-rate slip
x = zetamesh(creepindex);
x = (x - min(x));
cracksliprate = Vplate.*sqrt(1 - (x./w-1).^2);
sliprate(creepindex)=cracksliprate;
sliprate(zetamesh > (xlock+w)) = Vplate;
end