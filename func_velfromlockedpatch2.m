function dpred = func_velfromlockedpatch2(m,params)
% params should have the following: 
% Elastic parameters: mu,nu,
% station locations: ox,
% flags for which stations are [velx,velz]
%
% model parameters for inversion:
% m(1) - downdip locking limit in [km]
% m(2) - transition width (from lock to creep at plate rate) [km]
% m(3) - Vplate in [mm/yr]
% m(4) - fault dip in degrees
% m(5) - plate thickness in [km]
% AUTHOR:
% Rishav Mallick, JPL, 2024

% construct fault parameter vector for meshing
% [Vpl,x0,z0,W,δ,dw]
meshsize = 0.5e3;
fparameters = [1,0,1,300e3,m(4),meshsize;...
            -sind(m(4)),0,1,m(5)*1e3,m(4)+90,m(5)*1e3];
Vplate = m(3); % plate velocity in [mm/yr]
% station locations
ox = params.ox;
% create fault mesh
rcv = geometry.receiver(fparameters,geometry.LDhs(params.mu,params.nu));

% find patches on the megathrust (not the hinge) 
index = rcv.Vpl>0;
% compute distance along the megathrust in [km]
dipdist = sqrt(rcv.xc(:,1).^2 + rcv.xc(:,2).^2)./1e3;

% locked patches on the fault & the hinge too
Ilockedpatch = (dipdist<m(1) & index);

% enforce plate rate below some depth
Icreep = (dipdist > (m(1)+m(2)) | ~index);

% imposed BCs on locked & free-creep patches
IimposedBC = Ilockedpatch|Icreep;

sliprate=zeros(rcv.N,1);
% using the traction free boundaries for the stress-driven creeping
% segments while the forced slip segments are either no slip or full
% plate-rate slip
x = dipdist(~IimposedBC);
x = (x - min(x));
cracksliprate = Vplate.*sqrt(1 - (x./m(2)-1).^2);

sliprate(~IimposedBC)=cracksliprate;
sliprate(Icreep)=Vplate.*rcv.Vpl(Icreep);
sliprate(~index) = 0; % the hinge is not creeping but we don't want to deal with stresses from the hinge either

% construct appropriate displacement kernel for all stations
[Gx,Gz,~,~] = geometry.computeFaultDisplacementKernels(rcv,[ox,zeros(size(ox))]);
stationindex = (params.flag == 1);
bigG = [Gz(stationindex,:);Gx(~stationindex,:)];

% return predicted velocity field 
dpred = bigG*(sliprate-rcv.Vpl.*Vplate);
% for horizontals we need to account for reference frame and block motions
dpred(params.flag==0) = -dpred(params.flag==0) + Vplate.*(ox(params.flag==0)>0);
end
