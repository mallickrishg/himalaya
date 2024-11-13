function dpred = func_velfromlockedpatch(m,params)
% params should have the following: 
% Elastic parameters: mu,nu,
% station locations: ox,
% flags for which stations are [velx,velz]
%
% model parameters for inversion:
% m(1) - updip locking limit in [km]
% m(2) - downdip locking limit in [km]
% m(3) - transition width (from lock to creep at plate rate) [km]
% m(4) - Vplate in [mm/yr]
% m(5) - fault dip in degrees
% m(6) - plate thickness in [km]
% AUTHOR:
% Rishav Mallick, JPL, 2024

% construct fault parameter vector for meshing
% [Vpl,x0,z0,W,δ,dw]
meshsize = 0.5e3;
fparameters = [1,0,1,300e3,m(5),meshsize;...
            -sind(m(5)),0,1,m(6)*1e3,m(5)+90,m(6)*1e3];
Vplate = m(4); % plate velocity in [mm/yr]
% station locations
ox = params.ox;
% create fault mesh
rcv = geometry.receiver(fparameters,geometry.LDhs(params.mu,params.nu));

% find patches on the megathrust (not the hinge) 
index = rcv.Vpl>0;
% compute distance along the megathrust in [km]
dipdist = sqrt(rcv.xc(:,1).^2 + rcv.xc(:,2).^2)./1e3;

% locked patches on the fault & the hinge too
Ilockedpatch = (dipdist>m(1) & dipdist<m(2) & index);

% enforce plate rate below some depth
Icreep = (dipdist > (m(2)+m(3)) | ~index);

% imposed BCs on locked & free-creep patches
IimposedBC = Ilockedpatch|Icreep;

% compute traction and displacement kernels
[K,~,~,~] = geometry.computeFaultTractionKernels(rcv,rcv);
Vplvector = rcv.Vpl.*Vplate;
Vplvector(~index) = 0; % we don't want to account for stress from the hinge
% long-term stressing rate
sig0 = -K*(Vplvector);
[Gx,Gz,~,~] = geometry.computeFaultDisplacementKernels(rcv,[ox,zeros(size(ox))]);

% construct appropriate displacement kernel for all stations
stationindex = (params.flag == 1);
bigG = [Gz(stationindex,:);Gx(~stationindex,:)];

subK = K(~IimposedBC,~IimposedBC);
%alter the stressing rate by subtracting the stresses that develop from previously locking 
substressrate = sig0(~IimposedBC) - K(~IimposedBC,Icreep)*(Vplvector(Icreep));

subsliprate=subK\substressrate;
sliprate=zeros(rcv.N,1);
% using the traction free boundaries for the stress-driven creeping
% segments while the forced slip segments are either no slip or full
% plate-rate slip
sliprate(~IimposedBC)=subsliprate;
sliprate(Icreep)=Vplate.*rcv.Vpl(Icreep);
sliprate(~index) = 0; % the hinge is not creeping but we don't want to deal with stresses from the hinge either

% return predicted velocity field 
dpred = bigG*(sliprate-rcv.Vpl.*Vplate);
% for horizontals we need to account for reference frame and block motions
dpred(params.flag==0) = -dpred(params.flag==0) + Vplate.*(ox(params.flag==0)>0);
end
