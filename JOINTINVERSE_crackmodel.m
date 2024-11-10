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
Vplate = 18; % convergence on MHT in mm/yr
earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver(faultfilename,earthModel);

% load data [x(km),vz(mm/yr),σz(mm/yr)]
xpred = linspace(-50,300,500)'.*1e3;% predicted locations
[ox1,d1,Cd1,flag1] = create_inputdataset('data/InSAR_vel_profile.txt','vertical');
[ox2,d2,Cd2,flag2] = create_inputdataset('data/fpp_panda_narrow.dat','horizontal');

% convert horizontals into upper-plate reference frame
d2 = Vplate - d2;
d2(ox2<0) = d2(ox2<0) - Vplate;

% combine datasets and create data covariance matrix
reweight = 1.0;
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
% m(3) - down-dip free-slip extent in [km]
% m(4) - Vplate
% m(5) - fault dip
% m(6) - plate thickness in [km]

rng(42)
Nsamples = 1e2;

parameters = [];
parameters.ox = ox;
parameters.flag = flag;
parameters.mu = mu;
parameters.nu = nu;

% provide model initialization: ideally this is an optimization solution 
minit = [1,100,130,20,10,100];

% non-linear sampling with slice sampler
dpred = @(m) func_velfromlockedpatch(m,parameters);
residuals=@(dpred,d,W) sum(diag(W))/(sum(diag(W))^2 - sum(diag(W).^2))*(dpred-d)'*W*(dpred-d);
% setup priors (bounds)

LB =  [0,10,10,14,5,1];
UB =  [100,500,500,25,20,500];
mnames = {'\zeta_{up-dip}';'\zeta_{down-dip}';'\zeta_{free}';'V_{pl}';'\delta';'T_{plate}'};

% log-likelihood on a transformed variable (makes sampling easier)
likelihoodfun = @(m) -0.5.*(residuals(dpred((UB-LB)./(1+exp(-m)) + LB),d,W)) - sum(m + 2*log(1+exp(-m)));

tic
disp('Begin Slicesampling')

zetasamples = slicesample(minit,Nsamples,'logpdf',likelihoodfun,'thin',2,'burnin',Nsamples/100,'width',1);

mposterior = zeros(size(zetasamples));
for i = 1:length(LB)
    mposterior(:,i) = (UB(i)-LB(i))./(1+exp(-zetasamples(:,i))) + LB(i);
end

toc
disp('End sampling')

%% make posterior distribution figure
figure(1),clf
plot_joint_post_pdf(mposterior,mnames);





