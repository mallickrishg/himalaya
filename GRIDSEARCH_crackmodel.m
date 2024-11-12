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
[ox1,d1,Cd1,flag1] = create_inputdataset('data/InSAR_vel_profile.txt','vertical');
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
minit = [100,10,16,10,100];

% non-linear sampling with slice sampler
dpred = @(m) func_velfromlockedpatch2(m,parameters);
residuals=@(dpred,d,W) sum(diag(W))/(sum(diag(W))^2 - sum(diag(W).^2))*(dpred-d)'*W*(dpred-d);
% setup priors (bounds)

LB =  [50,5,15,8,20];
UB =  [200,100,25,15,200];
mnames = {'\zeta_{down-dip}';'W';'V_{pl}';'\delta';'T_{plate}'};

% compute optimization solution
options = optimoptions('lsqnonlin','Algorithm','interior-point',...
    'TypicalX',[5,5,0.5,1,5],'StepTolerance',1e-3,...
    'FiniteDifferenceStepSize',[1,1,1e-6,1e-6,1],'Display','none');
mopt = lsqnonlin(@(m) sqrt(W)*(dpred(m) - d),minit,LB,UB,[],[],[],[],[],options);
disp([minit' mopt'])

% compute predictions from optimization
params = [];
params.ox = xpred;
params.mu = mu;
params.nu = nu;
params.flag = zeros(length(xpred),1);
vxopt = func_velfromlockedpatch2(mopt,params);
params.flag = ones(length(xpred),1);
vzopt = func_velfromlockedpatch2(mopt,params);

%% execute gridsearch
% grid up domain
nglock = 20;
ngw = 20;
ngvpl = 20;
x1vec = linspace(80,160,nglock);
x2vec = linspace(0,100,ngw);
x3vec = linspace(15,25,ngvpl);
% N-D grid
[x1g,x2g,x3g] = meshgrid(x1vec,x2vec,x3vec);

misfit = zeros(size(x1g));
minmisfit = zeros(ngvpl,1);

% define loss function
dpred_gridsearch = @(m) func_velfromlockedpatch2([m,mopt(4:5)],parameters);
tic
for k = 1:length(x3vec)
    disp(['Completed ' num2str(round(k*100/length(x3vec))) ' %'])
    for i = 1:length(x2vec)
        parfor j = 1:length(x1vec)
            misfit(i,j,k) = residuals(dpred_gridsearch([x1g(i,j,k),x2g(i,j,k),x3g(i,j,k)]),d,W);
        end
    end
    dummy = misfit(:,:,k);
    minmisfit(k) = min(dummy(:));
end

toc
%% plot relative probability of models as marginal PDFs

figure(1),clf
dummy = exp(-(misfit-min(misfit(:))));

for i = 1:3    
    for j = 1:3
        if (i==j)
            subplot(3,3,(3*(i-1) + j))
            if i==1
                bar(x1vec,max(squeeze(max(dummy,[],1)),[],2),'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
                set(gca,'FontSize',15,'Color','none','XTick',[])
            elseif (i==2)
                bar(x2vec,max(squeeze(max(dummy,[],2)),[],2),'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
                set(gca,'FontSize',15,'Color','none','XTick',[])
            elseif (i==3)
                bar(x3vec,max(squeeze(max(dummy,[],1)),[],1),'FaceColor',rgb('skyblue'),'EdgeColor',rgb('skyblue'))
                axis tight, grid on
                xlabel('V_{pl} [mm/yr]')
                set(gca,'FontSize',15,'Color','none')
            else
                error('check number of variables to plot')
            end
            
        elseif i>j
            subplot(3,3,(3*(i-1) + j))
            if (3*(i-1) + j)==4
                pcolor(x1vec,x2vec,squeeze(max(dummy,[],3))), shading flat
                axis tight, box on, grid on
                alpha 0.6
                ylabel('W [km]')
                set(gca,'FontSize',15,'Color','none','XTick',[])
            elseif (3*(i-1) + j)==7
                pcolor(x1vec,x3vec,squeeze(max(dummy,[],1))'), shading flat
                axis tight, box on, grid on
                alpha 0.6
                xlabel('\zeta_{lock} [km]')
                ylabel('V_{pl} [mm/yr]')
                set(gca,'FontSize',15,'Color','none')
            elseif (3*(i-1) + j)==8
                pcolor(x2vec,x3vec,squeeze(max(dummy,[],2))'), shading flat
                axis tight, box on, grid on
                alpha 0.6
                xlabel('W [km]')
                set(gca,'FontSize',15,'Color','none','YTick',[])
            end
        end
        
    end
    
end
