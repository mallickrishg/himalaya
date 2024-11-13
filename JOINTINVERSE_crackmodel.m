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

% load data [x(km),vz(mm/yr),σz(mm/yr)]
xpred = linspace(-50,300,500)'.*1e3;% predicted locations
[ox1,d1,Cd1,flag1] = create_inputdataset('data/InSAR_vel_profile_2k.txt','vertical');
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
Nchains = 12;

parameters = [];
parameters.ox = ox;
parameters.flag = flag;
parameters.mu = mu;
parameters.nu = nu;

% provide model initialization: ideally this is an optimization solution 
minit = [100,10,16,10,100];

% non-linear sampling with slice sampler
dpred = @(m) func_velfromlockedpatch2(m,parameters);
% residuals=@(dpred,d,W) sum(diag(W))/(sum(diag(W))^2 - sum(diag(W).^2))*(dpred-d)'*W*(dpred-d);
residuals=@(dpred,d,W) (dpred-d)'*W*(dpred-d);

% setup priors (bounds)
LB =  [50,5,15,8,20];
UB =  [200,100,25,15,200];
mnames = {'\zeta_{lock} [km]';'W [km]';'V_{pl} [mm/yr]';'\delta º';'T_{plate} [km]'};

% compute optimization solution
options = optimoptions('lsqnonlin','Algorithm','interior-point',...
    'TypicalX',[5,5,0.5,1,5],'StepTolerance',1e-3,...
    'FiniteDifferenceStepSize',[1,1,1e-6,1e-6,1],'Display','none');
mopt = lsqnonlin(@(m) sqrt(W)*(dpred(m) - d),minit,LB,UB,[],[],[],[],[],options);
disp([minit' mopt'])

% log-likelihood on a logit-transformed variable (makes sampling easier)
likelihoodfun = @(m) -0.5.*(residuals(dpred((UB-LB)./(1+exp(-m)) + LB),d,W)) - sum(m + 2*log(1+exp(-m)));

tic
disp('Begin Slicesampling')

% sample with slicesample and parallelize over Nchains
zetasamples3d = zeros(Nsamples,Nchains,length(LB));
zetasamples = zeros(Nsamples*Nchains,length(LB));
mposterior = zeros(size(zetasamples));
zetainit = log(mopt-LB)-log(UB-mopt);
parfor chains = 1:Nchains
    zetasamples3d(:,chains,:) = slicesample(zetainit,Nsamples,...
        'logpdf',likelihoodfun,'thin',3,'burnin',Nsamples,'width',[3,1,1,1,1]);
end
% combine Nchains into a single matrix
for i = 1:Nchains
    zetasamples(Nsamples*(i-1) + 1:Nsamples*i,:) = squeeze(zetasamples3d(:,i,:));
end
% transform from logit coordinates to natural coordinates
for i = 1:length(LB)
    mposterior(:,i) = (UB(i)-LB(i))./(1+exp(-zetasamples(:,i))) + LB(i);
end

toc
disp('End sampling')

%% make posterior distribution figure
figure(1),clf
set(gcf,'Color','w')
plot_joint_post_pdf(mposterior,mnames);

% plot predictions
Npred = 301;
sampleid = randperm(Nchains*Nsamples,Npred);
vxpred = zeros(length(xpred),Npred);
vzpred = zeros(length(xpred),Npred);

parfor i = 1:length(sampleid)
    params = [];
    params.ox = xpred;
    params.flag = zeros(length(xpred),1);
    params.mu = mu;
    params.nu = nu;
    vxpred(:,i) = func_velfromlockedpatch2(mposterior(sampleid(i),:),params);
    params.flag = ones(length(xpred),1);
    vzpred(:,i) = func_velfromlockedpatch2(mposterior(sampleid(i),:),params);
end
% compute predictions from optimization
params = [];
params.ox = xpred;
params.mu = mu;
params.nu = nu;
params.flag = zeros(length(xpred),1);
vxopt = func_velfromlockedpatch2(mopt,params);
params.flag = ones(length(xpred),1);
vzopt = func_velfromlockedpatch2(mopt,params);

nmesh = 1000;
xmesh = zeros(nmesh,Npred);
slipratepred = zeros(nmesh,Npred);
parfor i = 1:length(sampleid)
    [xm,vm] = func_slipratefromlockedpatch(mposterior(sampleid(i),:),nmesh);
    xmesh(:,i) = xm;
    slipratepred(:,i) = vm;
end

% compute potency accumulation rate (this is per unit length along-strike)
potencyrate = zeros(Npred,1);
for i = 1:Npred
    potencyrate(i,1) = -trapz(xmesh(:,i)./cosd(mposterior(sampleid(i),4)),...
        slipratepred(:,i)-mposterior(sampleid(i),3)');
end
%% plot data predictions and slip rate on fault
figure(2),clf
set(gcf,'Color','w')
subplot(3,1,1)
errorbar(ox1./1e3,d1,sqrt(diag(Cd1)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vzpred,'-','Color',[1 0 0 0.05])
plot(xpred./1e3,vzopt,'k-','Linewidth',2)
axis tight
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(3,1,2)
errorbar(ox2./1e3,d2,sqrt(diag(Cd2)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vxpred,'-','Color',[1 0 0 0.05])
plot(xpred./1e3,vxopt,'k-','Linewidth',2)
axis tight
xlim([-50,300])
xlabel('x (km)'), ylabel('v_x [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

% plot fault slip rate
subplot(3,1,3)
plot(xmesh./1e3,slipratepred,'-','Color',[1 0 0 0.1]), hold on
[xm,vm] = func_slipratefromlockedpatch(mopt,nmesh);
plot(xm./1e3,vm,'k-','Linewidth',2)
xlabel('x (km)'), ylabel('sliprate [mm/yr]')
xlim([-50,300])
ylim([0 25])
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

% potency accumulation rate
figure(3),clf
histogram(potencyrate,'Normalization','probability')
axis tight
xlabel('potency rate [mm/yr-m]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both','XScale','log')
