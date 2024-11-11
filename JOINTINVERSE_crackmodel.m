% Script to read vertical velocities and carry out an unregularized inverse
% problem using a crack-based interseismic model following 
% Mallick et al., (2021)
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear
addpath functions/
import('geometry.*')

% Elastic parameters (homogenous medium)
nu = 0.25;% Poisson's ratio
mu = 30e3;% in MPa
faultfilename = 'megathrust2d.seg';

% load megathrust mesh
earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver(faultfilename,earthModel);

% load data [x(km),vz(mm/yr),σz(mm/yr)]
xpred = linspace(-50,300,500)'.*1e3;% predicted locations
[ox1,d1,Cd1,flag1] = create_inputdataset('data/InSAR_vel_profile.txt','vertical');
[ox2,d2,Cd2,flag2] = create_inputdataset('data/fpp_panda.dat','horizontal');

% combine datasets and create data covariance matrix
ox = [ox1;ox2];
d = [d1;d2];
Cd = blkdiag(Cd1.*length(d)/length(d1),Cd2.*length(d)/length(d2));
flag = [flag1;flag2];

% weighting function is inverse of data covariance matrix
W = inv(Cd);
%% construct loss function for inversion
% model parameters for inversion:
% m(1) - updip locking limit in [km]
% m(2) - downdip locking limit in [km]
% m(3) - transition width (from lock to creep at plate rate) [km]
% m(4) - Vplate in [mm/yr]
% m(5) - fault dip in degrees
% m(6) - plate thickness in [km]

rng(42)
Nsamples = 1e2;

parameters = [];
parameters.ox = ox;
parameters.flag = flag;
parameters.mu = mu;
parameters.nu = nu;

% provide model initialization: ideally this is an optimization solution 
minit = [10,100,10,20,10,100];

% non-linear sampling with slice sampler
dpred = @(m) func_velfromlockedpatch(m,parameters);
residuals=@(dpred,d,W) sum(diag(W))/(sum(diag(W))^2 - sum(diag(W).^2))*(dpred-d)'*W*(dpred-d);
% setup priors (bounds)

LB =  [0,50,5,15,8,90];
UB =  [50,200,100,25,15,110];
mnames = {'\zeta_{up-dip}';'\zeta_{down-dip}';'W';'V_{pl}';'\delta';'T_{plate}'};

% compute optimization solution
mopt = lsqnonlin(@(m) sqrt(W)*(dpred(m) - d),minit,LB,UB);
disp(mopt)

% log-likelihood on a transformed variable (makes sampling easier)
likelihoodfun = @(m) -0.5.*(residuals(dpred((UB-LB)./(1+exp(-m)) + LB),d,W)) - sum(m + 2*log(1+exp(-m)));

tic
disp('Begin Slicesampling')

zetasamples = slicesample(mopt,Nsamples,'logpdf',likelihoodfun,'thin',2,'burnin',Nsamples/10,'width',[3,3,3,3,2,3]);

mposterior = zeros(size(zetasamples));
for i = 1:length(LB)
    mposterior(:,i) = (UB(i)-LB(i))./(1+exp(-zetasamples(:,i))) + LB(i);
end

toc
disp('End sampling')

%% make posterior distribution figure
figure(1),clf
plot_joint_post_pdf(mposterior,mnames);

% plot predictions
vxpred = zeros(length(xpred),Nsamples);
vzpred = zeros(length(xpred),Nsamples);
for i = 1:10:Nsamples
    params = [];
    params.ox = xpred;
    params.flag = zeros(length(xpred),1);
    params.mu = mu;
    params.nu = nu;
    vxpred(:,i) = func_velfromlockedpatch(mposterior(i,:),params);
    params.flag = ones(length(xpred),1);
    vzpred(:,i) = func_velfromlockedpatch(mposterior(i,:),params);
end
% compute predictions from optimization
params.flag = zeros(length(xpred),1);
vxopt = func_velfromlockedpatch(mopt,params);
params.flag = ones(length(xpred),1);
vzopt = func_velfromlockedpatch(mopt,params);

%% plot data predictions and fault slip rate
figure(2),clf
set(gcf,'Color','w')
subplot(2,1,1)
errorbar(ox1./1e3,d1,sqrt(diag(Cd1)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vzpred,'-','Color',[1 0 0 0.2])
plot(xpred./1e3,vzopt,'r-','Linewidth',2)
axis tight
% ylim([-1,1]*5)
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(2,1,2)
errorbar(ox2./1e3,d2,sqrt(diag(Cd2)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vxpred,'-','Color',[1 0 0 0.1])
plot(xpred./1e3,vxopt,'r-','Linewidth',2)
axis tight
% ylim([-0.2,1.2])
xlim([-50,300])
xlabel('x (km)'), ylabel('v_x/v_{plate}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

% subplot(3,1,3)
% plot(rcv.xc(:,1)./1e3,msol + rcv.Vpl*Vplate,'ks','LineWidth',2), hold on
% scatter(rcv.xc(:,1)./1e3,msamples + rcv.Vpl*Vplate,20,[1,0,0],'filled'), alpha 0.1
% hold on
% plot(rcv.xc(:,1)./1e3,median(msamples,2) + rcv.Vpl*Vplate,'ro','LineWidth',2,'MarkerFaceColor','r')
% xlabel('x (km)'), ylabel('slip rate [mm/yr]')
% axis tight, box on
% xlim([-50,300])
% set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

%% testing
minit = [20,30,50,20,10,100];
options = optimoptions('lsqnonlin','Algorithm','interior-point',...
    'TypicalX',[5,10,5,0.5,1,5],'StepTolerance',2);
% compute optimization solution
mopt = lsqnonlin(@(m) sqrt(W)*(dpred(m) - d),minit,LB,UB);
disp(mopt')
% compute predictions from optimization
params.flag = zeros(length(xpred),1);
vxopt = func_velfromlockedpatch(mopt,params);
params.flag = ones(length(xpred),1);
vzopt = func_velfromlockedpatch(mopt,params);

figure(2),clf
set(gcf,'Color','w')
subplot(2,1,1)
errorbar(ox1./1e3,d1,sqrt(diag(Cd1)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vzopt,'r-','Linewidth',2)
axis tight
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(2,1,2)
errorbar(ox2./1e3,d2,sqrt(diag(Cd2)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vxopt,'r-','Linewidth',2)
axis tight
xlim([-50,300])
xlabel('x (km)'), ylabel('v_x/v_{plate}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')