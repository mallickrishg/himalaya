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

% plot geometry
figure(2),clf
plotpatch2d(rcv), hold on
quiver(rcv.xc(:,1)./1e3,rcv.xc(:,2)./1e3,rcv.nv(:,1),rcv.nv(:,2),'r')
quiver(rcv.xc(:,1)./1e3,rcv.xc(:,2)./1e3,rcv.dv(:,1),rcv.dv(:,2),'b')
axis tight equal
colorbar
clim([-1 1])
colormap("turbo(10)")