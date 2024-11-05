% script to visualize surface velocity from ESPM (Kanda & Simons 2010)
% AUTHOR: 
% Rishav Mallick, JPL, 2024

clear
addpath functions/
import('geometry.*')

% provide distance from trench in [km]
x = linspace(-100,300,1e3)';
% provide reference point in [km] (large x value +ve value indicates upper plate)
x0 = 1e6;

% m - vector containing [Vpl,dip(º),locking depth(km),plate thickness(km)]
m = [1,8,20,50];

% compute velocity field
[vx,vz] = espm2d(m,x);
% compute hoirzontal velocity at reference point
[vx0,~] = espm2d(m,x0);

figure(1),clf
subplot(2,1,1)
plot(x,-(vx-vx0),'LineWidth',2)
axis tight
ylim([0 1])
xlabel('x (km)'), ylabel('v_x/v_{pl}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(2,1,2)
plot(x,vz,'LineWidth',2)
axis tight
ylim([-1 1]*max(abs(vz)))
xlabel('x (km)'), ylabel('v_z/v_{pl}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

%% read mesh file and plot internal displacements
% Elastic parameters (homogenous medium)
nu=0.25;% Poisson's ratio
mu=1;% in MPa

% load megathrust mesh
earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver('megathrust2d.seg',earthModel);

% provide observation points
nx = 200;
nz = 201;
xobs = linspace(-10,301,nx).*1e3;
zobs = linspace(-60,0,nz).*1e3;
[xg,zg] = meshgrid(xobs,zobs);
ox = [xg(:),zg(:)];
% compute displacement kernels
[Gdx,Gdz,~,~] = geometry.computeFaultDisplacementKernels(rcv,ox);

% provide a slip distrubtion for the mesh
slip = -rcv.Vpl; %slip(rcv.Vpl == 1) = 0;
ux = Gdx*slip;
uz = Gdz*slip;

% plot geometry
nskip = 19;
figure(2),clf
subplot(3,1,[1 2])
plotpatch2d(rcv), hold on, axis tight equal
% quiver(rcv.xc(:,1)./1e3,rcv.xc(:,2)./1e3,rcv.nv(:,1),rcv.nv(:,2),'r')
% quiver(rcv.xc(:,1)./1e3,rcv.xc(:,2)./1e3,rcv.dv(:,1).*rcv.Vpl,rcv.dv(:,2).*rcv.Vpl,'b')
% scatter(xobs(:)./1e3,zobs(:)./1e3,50,uz,'filled')
pcolor(xobs./1e3,zobs./1e3,reshape(uz,nz,nx)), shading interp
quiver(xg(1:nskip:end)'./1e3,zg(1:nskip:end)'./1e3,ux(1:nskip:end),uz(1:nskip:end),'k')
axis tight equal
cb=colorbar;cb.Label.String='v_z';cb.Location='northoutside';
clim([-1 1]*0.25)
colormap("bluewhitered(100)")
xlabel('x (km)'), ylabel('depth (km)')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(3,1,3)
index = zg(:) == 0;
plot(xg(index)./1e3,uz(index),'o-'), hold on
plot(xg(index)./1e3,ux(index),'o-')
axis tight
xlabel('x (km)'), ylabel('v_z/v_{pl}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')