% Script to read vertical velocities and carry out an unregularized inverse
% problem accounting for data and model uncertainties following: 
% Duputel et al., 2014 & Ragon et al., 2019
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear
addpath functions/
import('geometry.*')

% Elastic parameters (homogenous medium)
nu = 0.25;% Poisson's ratio
mu = 1;% in MPa
faultfilename = 'megathrust2d.seg';

% load megathrust mesh
Vplate = 20; % convergence on MHT in mm/yr
earthModel = geometry.LDhs(mu,nu);
rcv = geometry.receiver(faultfilename,earthModel);

% load data [x(km),vz(mm/yr),σz(mm/yr)]
xpred = linspace(-50,300,500)'.*1e3;% predicted locations
[ox1,d1,Cd1,flag1] = create_inputdataset('data/InSAR_vel_profile.txt','vertical');
[ox2,d2,Cd2,flag2] = create_inputdataset('data/fpp_panda.dat','horizontal');

% convert horizontals into upper-plate reference frame
d2 = Vplate - d2;
d2(ox2<0) = d2(ox2<0) - Vplate;

% combine datasets and create data covariance matrix
reweight = 1.5;
ox = [ox1;ox2];
d = [d1;d2];
Cd = blkdiag(Cd1.*reweight,Cd2.*(2-reweight));
flag = [flag1;flag2];

%% compute Greens functions and assemble design matrix for inverse problem
[Gx,Gz,~,~] = geometry.computeFaultDisplacementKernels(rcv,[ox,zeros(size(ox))]);
[Gxpred,Gzpred,~,~] = geometry.computeFaultDisplacementKernels(rcv,[xpred,zeros(size(xpred))]);

index = (flag == 1);
bigG = [Gz(index,:);Gx(~index,:)];

%% compute bounded-value least squares optimization
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
msol = lsqlin(sqrtm(inv(Cd))*bigG,sqrtm(inv(Cd))*d,[],[],[],[],lb,ub);

%% sample from posterior distribution using HMC sampling
Nsamples = 100;
alpha2 = 1e-6; % this is a relative weight parameter that just needs to be small enough to promote sampling speed, but not affect results
M = bigG'*(Cd\bigG) + alpha2*eye(rcv.N);
r = (d'*(Cd\bigG))';
F = [eye(rcv.N);-eye(rcv.N)];
g = [-lb;ub] + 1e-9;

% compute posterior distribution
[msamples, ~] = HMC_exact(F, g, M, r, 0, Nsamples, msol);

% compute residual pdf
residuals = d - bigG*msamples;
residuals_pdf = zeros(Nsamples,1);
for i = 1:Nsamples
    residuals_pdf(i) = residuals(:,i)'*(Cd\residuals(:,i));
end

%% plot residual PDF
figure(10),clf
histogram(sqrt(residuals_pdf./length(d)),'Normalization','pdf','EdgeColor','none'), hold on
plot(mean(sqrt(diag(Cd1))).*[1 1],get(gca,'YLim'),'k-','LineWidth',2)
plot(mean(sqrt(diag(Cd2))).*[1 1],get(gca,'YLim'),'r-','LineWidth',2)
legend('pdf','\sigma_{insar}','\sigma_{gnss}','location','best','box','off')
xlabel('\sigma [mm/yr]'), ylabel('pdf')
xlim([0 2.5])
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')
%% plot results
figure(1),clf
set(gcf,'Color','w')
subplot(3,1,1)
errorbar(ox1./1e3,d1,sqrt(diag(Cd1)),'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,Gzpred*msol,'k-','LineWidth',2)
plot(xpred./1e3,Gzpred*msamples(:,2:end),'-','Color',[1 0 0 0.1])
axis tight
ylim([-1,1]*5)
xlabel('x (km)'), ylabel('v_z [mm/yr]')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(3,1,2)
errorbar(ox2./1e3,(d2 + 0.*Vplate.*(ox2<0))./Vplate,sqrt(diag(Cd2))./Vplate,'o','LineWidth',1,'CapSize',0,'MarkerFaceColor','blue'), hold on
plot(xpred./1e3,(Gxpred*msol + 0.*Vplate.*(xpred<0))./Vplate,'k-','LineWidth',2)
plot(xpred./1e3,(Gxpred*msamples(:,2:end) + 0.*Vplate.*(xpred<0))./Vplate,'-','Color',[1 0 0 0.1])
axis tight
ylim([-0.2,1.2])
xlabel('x (km)'), ylabel('v_x/v_{plate}')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(3,1,3)
plot(rcv.xc(:,1)./1e3,msol + rcv.Vpl*Vplate,'k.','LineWidth',2), hold on
scatter(rcv.xc(:,1)./1e3,msamples + rcv.Vpl*Vplate,20,[1,0,0],'filled'), alpha 0.1
hold on
plot(rcv.xc(:,1)./1e3,median(msamples,2) + rcv.Vpl*Vplate,'ro','LineWidth',2,'MarkerFaceColor','r')
xlabel('x (km)'), ylabel('slip rate [mm/yr]')
axis tight, box on
xlim([-50,300])
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')
