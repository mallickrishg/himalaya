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

%% conduct gridsearch over select parameters
% model parameters for inversion:
% m(1) - downdip locking limit in [km]
% m(2) - transition width (from lock to creep at plate rate) [km]
% m(3) - Vplate in [mm/yr]
% m(4) - fault dip in degrees
% m(5) - plate thickness in [km]

parameters = [];
parameters.ox = ox;
parameters.flag = flag;
parameters.mu = mu;
parameters.nu = nu;

% provide model initialization: ideally this is an optimization solution 
minit = [10,100,10,16,10,100];

% non-linear sampling with slice sampler
dpred = @(m) func_velfromlockedpatch(m,parameters);
% residuals=@(dpred,d,W) sum(diag(W))/(sum(diag(W))^2 - sum(diag(W).^2))*(dpred-d)'*W*(dpred-d);
residuals=@(dpred,d,W) (dpred-d)'*W*(dpred-d);

% setup priors (bounds)
LB =  [0,80,5,15,5,20];
UB =  [80,200,100,25,15,200];
mnames = {'\zeta_{creep} [km]';'\zeta_{lock} [km]';'W [km]';'V_{pl} [mm/yr]';'\delta º';'T_{plate} [km]'};

% compute optimization solution
options = optimoptions('lsqnonlin','Algorithm','interior-point',...
    'TypicalX',[5,5,5,0.5,1,5],'StepTolerance',1e-3,...
    'FiniteDifferenceStepSize',[1,1,1,1e-6,1e-6,1],'Display','none');
mopt = lsqnonlin(@(m) sqrt(W)*(dpred(m) - d),minit,LB,UB,[],[],[],[],[],options);
disp([minit' mopt'])

% compute predictions from optimization
params = [];
params.ox = xpred;
params.mu = mu;
params.nu = nu;
params.flag = zeros(length(xpred),1);
vxopt = func_velfromlockedpatch(mopt,params);
params.flag = ones(length(xpred),1);
vzopt = func_velfromlockedpatch(mopt,params);

%% execute gridsearch
% grid up domain
ngup = 10;
nglock = 20;
ngw = 20;
ngvpl = 10;
ngdip = 10;
ngTe = 10;

x1vec = linspace(0,50,ngup);
x2vec = linspace(90,110,nglock);
x3vec = linspace(0,60,ngw);
x4vec = linspace(17,21,ngvpl);
x5vec = linspace(6,11,ngdip);
x6vec = linspace(30,100,ngTe);
% N-D grid
[x1g,x2g,x3g,x4g,x5g,x6g] = ndgrid(x1vec,x2vec,x3vec,x4vec,x5vec,x6vec);

misfit = zeros(size(x1g));

% define loss function
dpred_gridsearch = @(m) func_velfromlockedpatch(m,parameters);
tic
parfor k = 1:length(x4g(:))   
    misfit(k) = residuals(dpred_gridsearch([x1g(k),x2g(k),x3g(k),x4g(k),x5g(k),x6g(k)]),d,W);
end

toc
%% plot relative probability of models as marginal PDFs
figure(1),clf
set(gcf,'Color','w')
dummy = exp(-(misfit-min(misfit(:))));

for i = 1:length(size(misfit))
    for j = 1:length(size(misfit))
        if i==j
            subplot(6,6,(6*(i-1) + j))
            if i == 1
                toplot = squeeze(max(dummy,[],[6,5,4,3,2]));
                bar(x1vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
            elseif i == 2
                toplot = squeeze(max(dummy,[],[1,3,4,5,6]));
                bar(x2vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
            elseif i == 3
                toplot = squeeze(max(dummy,[],[1,2,4,5,6]));
                bar(x3vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
            elseif i == 4
                toplot = squeeze(max(dummy,[],[1,2,3,5,6]));
                bar(x4vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on 
            elseif i == 5
                toplot = squeeze(max(dummy,[],[1,2,3,4,6]));
                bar(x5vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on 
            elseif i == 6
                toplot = squeeze(max(dummy,[],[1,2,3,4,5]));
                bar(x6vec,toplot,'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
                xlabel(mnames{6})
            end
            set(gca,'FontSize',15,'Color','none')

        elseif i>j
            subplot(6,6,(6*(i-1) + j))
            if (6*(i-1) + j)==7
                toplot = squeeze(max(dummy,[],[6,5,4,3]))';
                pcolor(x1vec,x2vec,toplot), shading interp
                ylabel(mnames{2}) 

            elseif (6*(i-1) + j)==13
                toplot = squeeze(max(dummy,[],[6,5,4,2]))';
                pcolor(x1vec,x3vec,toplot), shading interp
                ylabel(mnames{3})
            elseif (6*(i-1) + j)==14
                toplot = squeeze(max(dummy,[],[6,5,4,1]))';
                pcolor(x2vec,x3vec,toplot), shading interp

            elseif (6*(i-1) + j)==19
                toplot = squeeze(max(dummy,[],[6,5,3,2]))';
                pcolor(x1vec,x4vec,toplot), shading interp
                ylabel(mnames{4})
            elseif (6*(i-1) + j)==20
                toplot = squeeze(max(dummy,[],[6,5,3,1]))';
                pcolor(x2vec,x4vec,toplot), shading interp
            elseif (6*(i-1) + j)==21
                toplot = squeeze(max(dummy,[],[6,5,1,2]))';
                pcolor(x3vec,x4vec,toplot), shading interp

            elseif (6*(i-1) + j)==25
                toplot = squeeze(max(dummy,[],[6,4,3,2]))';
                pcolor(x1vec,x5vec,toplot), shading interp
                ylabel(mnames{5})
            elseif (6*(i-1) + j)==26
                toplot = squeeze(max(dummy,[],[6,4,3,1]))';
                pcolor(x2vec,x5vec,toplot), shading interp
            elseif (6*(i-1) + j)==27
                toplot = squeeze(max(dummy,[],[6,4,1,2]))';
                pcolor(x3vec,x5vec,toplot), shading interp
            elseif (6*(i-1) + j)==28
                toplot = squeeze(max(dummy,[],[6,3,2,1]))';
                pcolor(x4vec,x5vec,toplot), shading interp

            elseif (6*(i-1) + j)==31
                toplot = squeeze(max(dummy,[],[5,4,3,2]))';
                pcolor(x1vec,x6vec,toplot), shading interp
                xlabel(mnames{1})
                ylabel(mnames{6})
            elseif (6*(i-1) + j)==32
                toplot = squeeze(max(dummy,[],[5,4,3,1]))';
                pcolor(x2vec,x6vec,toplot), shading interp
                xlabel(mnames{2})
            elseif (6*(i-1) + j)==33
                toplot = squeeze(max(dummy,[],[5,4,1,2]))';
                pcolor(x3vec,x6vec,toplot), shading interp
                xlabel(mnames{3})
            elseif (6*(i-1) + j)==34
                toplot = squeeze(max(dummy,[],[5,3,2,1]))';
                pcolor(x4vec,x6vec,toplot), shading interp
                xlabel(mnames{4})
            elseif (6*(i-1) + j)==35
                toplot = squeeze(max(dummy,[],[4,3,2,1]))';
                pcolor(x5vec,x6vec,toplot), shading interp
                xlabel(mnames{5})
            end
            alpha 0.6
            axis tight, box on, grid on
            set(gca,'FontSize',15,'Color','none')
        end
    end
end
colormap(turbo(20))
% return
%% plot results
% mopt(end) = 50;
nmesh = 300;

in = dummy>=0.5;
x1plot = x1g(in);
x2plot = x2g(in);
x3plot = x3g(in);
x4plot = x4g(in);
x5plot = x5g(in);
x6plot = x6g(in);
vxpred = zeros(length(xpred),length(x1plot));
vzpred = zeros(length(xpred),length(x1plot));
slipratepred = zeros(nmesh,length(x1plot));
xmesh = zeros(nmesh,length(x1plot));

parfor i = 1:length(x1plot)
    params = [];
    params.ox = xpred;
    params.flag = zeros(length(xpred),1);
    params.mu = mu;
    params.nu = nu;
    modelparams = [x1plot(i),x2plot(i),x3plot(i),x4plot(i),x5plot(i),x6plot(i)];
    vxpred(:,i) = func_velfromlockedpatch(modelparams,params);
    params.flag = ones(length(xpred),1);
    vzpred(:,i) = func_velfromlockedpatch(modelparams,params);
    
    [xm,vm] = func_slipratefromlockedpatch(modelparams,nmesh,params);
    xmesh(:,i) = xm(1:end-1);
    slipratepred(:,i) = vm(1:end-1);
end

%% plot data predictions and fault slip rate
figure(2),clf
set(gcf,'Color','w')
subplot(3,1,1)
errorbar(ox1./1e3,d1,sqrt(diag(Cd1)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vzpred,'-','Color',[1 0 0 0.2])
% plot(xpred./1e3,vzopt,'k-','Linewidth',2)
axis tight
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(3,1,2)
errorbar(ox2./1e3,d2,sqrt(diag(Cd2)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,vxpred,'-','Color',[1 0 0 0.2])
% plot(xpred./1e3,vxopt,'k-','Linewidth',2)
axis tight
xlim([-50,300])
xlabel('x (km)'), ylabel('v_x [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

% plot fault slip rate
subplot(3,1,3)
plot(xmesh./1e3,slipratepred,'-','Color',[1 0 0 0.1]), hold on
xlabel('x (km)'), ylabel('sliprate [mm/yr]')
xlim([-50,300])
ylim([0 25])
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')