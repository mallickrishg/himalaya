clear
addpath functions/
import('geometry.*')

% Elastic parameters (homogenous medium)
nu = 0.25;% Poisson's ratio
mu = 1;% in MPa

% load megathrust mesh
Vplate = 20; % convergence on MHT in mm/yr
earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver('megathrust2d.seg',earthModel);

% load data [x(km),vz(mm/yr),σz(mm/yr)]
datainput = readmatrix('data/InSAR_vel_profile.txt');
ox = datainput(:,1).*1e3; % convert to [m]
xpred = linspace(-50,200,500)'.*1e3;% predicted locations

d = datainput(:,2);
Cd = diag(datainput(:,3).^2)./4;

%% construct design matrix and weighting function
% compute displacement kernels
[~,Gz,~,~] = geometry.computeFaultDisplacementKernels(rcv,[ox,zeros(size(ox))]);
[Gxpred,Gzpred,~,~] = geometry.computeFaultDisplacementKernels(rcv,[xpred,zeros(size(xpred))]);
% compute traction kernels for regularization
[Kdd,~,~,~] = geometry.computeFaultTractionKernels(rcv,rcv);

%% construct linear inverse problem and solve using regularization
alpha2 = 1e-4;
% m0 = -rcv.Vpl.*15;
% msol = (Gz'*(Cd\Gz) + (alpha^2).*(Kdd'*Kdd))\(Gz'*(Cd\d));
% msol = (Gz'*(Cd\Gz) + (alpha^2)*eye(rcv.N))\(Gz'*(Cd\d) + (alpha^2)*eye(rcv.N)*m0);

% bounded value optimization
lb = zeros(rcv.N,1);
ub = zeros(rcv.N,1);
% for megathrust backslip -Vpl < v < 0
index = rcv.Vpl>0; 
% for hinge backslip 0 < v < Vpl*sin(δ)
lb(index) = -rcv.Vpl(index)*Vplate;
lb(~index) = -rcv.Vpl(~index)*Vplate*0.9;
ub(~index) = -rcv.Vpl(~index)*Vplate;
% optimization
msol = lsqlin(sqrt(inv(Cd))*Gz,sqrt(inv(Cd))*d,[],[],[],[],lb,ub);

% solve using HMC sampling
Nsamples = 100;
M = Gz'*(Cd\Gz) + alpha2*eye(rcv.N);
r = (d'*(Cd\Gz))';
F = [eye(rcv.N);-eye(rcv.N)];
g = [-lb;ub] + 1e-12;

% compute posterior distribution
[msamples, ~] = HMC_exact(F, g, M, r, 0, Nsamples, msol);

%% plot results
figure(1),clf
errorbar(ox./1e3,d,datainput(:,3),'o','LineWidth',1), hold on
plot(xpred./1e3,Gzpred*msol,'k-','LineWidth',2)
plot(xpred./1e3,Gzpred*msamples(:,2:end),'-','Color',[1 0 0 0.2])
axis tight
ylim([-1,1]*5)
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

figure(2),clf
plot(rcv.xc(:,1)./1e3,msol + rcv.Vpl*Vplate,'o','LineWidth',2), hold on
scatter(rcv.xc(:,1)./1e3,msamples + rcv.Vpl*Vplate,20,[1,0,0],'filled'), alpha 0.1
plot(rcv.xc(:,1)./1e3,median(msamples,2) + rcv.Vpl*Vplate,'ro','LineWidth',2,'MarkerFaceColor','r')
xlabel('x (km)'), ylabel('slip rate [mm/yr]')
axis tight
xlim([-50,200])
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')