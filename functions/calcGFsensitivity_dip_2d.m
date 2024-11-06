function [Kx_dip,Kz_dip] = calcGFsensitivity_dip_2d(rcv,patchfname,ox,oz,dipdelta)
% % Calculate sensitivity of Greens functions from epistemic uncertainties (in hidden parameters) % %
% calcGFsensitivity_dip_2d will calculate a sensitivity kernel for a perturbation to model-space 
% (not within the linear problem) and its resulting effect in the observational domain.
%
% OUTPUT
% K is a nobs x Npatch matrix
% INPUTS
% rcv - unicycle fault object
% ox,oz - station coordinates
% dipdelta - perturbation to dip of rcv (provide a vector that is Npatches)
% References : Ragon et al., (2019) in GJI
% Rishav Mallick, Caltech Seismolab 2023

neval = 10;

G = @(x) [ones(size(x)) x];
bigGx = zeros(length(ox),neval);
bigGz = zeros(length(ox),neval);
Kx_dip = zeros(length(ox),neval);
Kz_dip = zeros(length(ox),neval);

for i = 1:rcv.N
    rcvi = geometry.receiver(patchfname,rcv.earthModel);
    dipvec = linspace(rcv.dip(i)-dipdelta(i),rcv.dip(i)+dipdelta(i),neval)';
    for count = 1:neval
        rcvi.dip(i) = dipvec(count);
        [Gdx,Gdz,~,~] = geometry.computeFaultDisplacementKernels(rcvi,[ox oz]);
        bigGx(:,count) = Gdx(:,i);
        bigGz(:,count) = Gdz(:,i);
    end
    
    mfitx = zeros(length(ox),2);
    mfitz = zeros(length(ox),2);
    for count = 1:length(ox)
        mfitx(count,:) = G(dipvec)\bigGx(count,:)';
        mfitz(count,:) = G(dipvec)\bigGz(count,:)';
    end
    
    Kx_dip(:,i) = mfitx(:,2);
    Kz_dip(:,i) = mfitz(:,2);
end


end