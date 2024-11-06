function [Kx_pos,Kz_pos] = calcGFsensitivity_position_2d(rcv,patchfname,ox,oz,xdelta)
% % Calculate sensitivity of Greens functions from epistemic uncertainties (in hidden parameters) % %
% calcGFsensitivity_position_2d will calculate a sensitivity kernel for a perturbation to model-space 
% (not within the linear problem) and its resulting effect in the observational domain.
%
% OUTPUT
% K is a nobs x Npatch matrix
% INPUTS
% rcv - unicycle fault object
% ox,oz - station coordinates
% xdelta - perturbation to position of rcv (provide a vector that is Npatches)
% References : Ragon et al., (2019) in GJI
% Rishav Mallick, Caltech Seismolab 2023

neval = 10;

G = @(x) [ones(size(x)) x];
bigGx = zeros(length(ox),neval);
bigGz = zeros(length(ox),neval);
Kx_pos = zeros(length(ox),neval);
Kz_pos = zeros(length(ox),neval);

for i = 1:rcv.N
    rcvi = geometry.receiver(patchfname,rcv.earthModel);
    posvec = linspace(-xdelta(i),xdelta(i),neval)';
    for count = 1:neval
        rcvi.xc(i,:) = rcv.xc(i,:) + posvec(count).*rcv.nv(i,:);
        [Gdx,Gdz,~,~] = geometry.computeFaultDisplacementKernels(rcvi,[ox oz]);
        bigGx(:,count) = Gdx(:,i);
        bigGz(:,count) = Gdz(:,i);
    end
    
    mfitx = zeros(length(ox),2);
    mfitz = zeros(length(ox),2);
    for count = 1:length(ox)
        mfitx(count,:) = G(posvec)\bigGx(count,:)';
        mfitz(count,:) = G(posvec)\bigGz(count,:)';
    end
    
    Kx_pos(:,i) = mfitx(:,2);
    Kz_pos(:,i) = mfitz(:,2);
end


end