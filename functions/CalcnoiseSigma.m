function [nrmsR, tauf, dchi] = CalcnoiseSigma(Data,res)
% Function to compute realistic sigmas
% [nrmsR, tauf, dchi] = CalcnoiseSigma(Data,res)
% Algorithm here uses the changes in the nrms with increasing averaging
% interval to predict infinite averaging NRMS
% Data [3xnobs]
% Data(1,:) = Time in decyr
% Data(2,:) = Position (use same units for all input values)
% Data(3,:) = Uncertainties/weights on position
% res = residual positions after fitting a linear combination of signals

% First get some bounds on how much averaging we can do.
times = Data(1,:); err = Data(3,:);
start = min(times); stop = max(times);
% calculate minmax of averaging window in days
minav = 7; maxav = (stop-start)*365.25/10;

if maxav < minav 
	minav = maxav/4 
end

numav = fix(maxav/minav);

% zero the arrays we need
savsta = []; timsta = [];
% Get the "white-noise" part of the error budget by differincing adjacent
% data points.
nd = length(res);
dres = res(1:nd-1)-res(2:nd);
dvar = err(1:nd-1).^2+err(2:nd).^2;
dchi = sqrt(sum(dres.^2./dvar)/nd);
fprintf('RealSigma white dchi %10.6f\n',dchi);
err = err*dchi;
% Loop over the averaging intervals
for i = 1:numav
    dyr = i*minav/365.25;
    % Compute the averaged residals for this interval
    chi2 = 0; num = 0;
    avR = []; avS = [];
    for t = start:dyr:stop
        sel = find(times >= t & times < t+dyr);
        ars = res(sel); ass = err(sel);
        if length(ars) > 2
            w = 1./ass.^2;
            Wmean = w*ars'/sum(w);
            Werr = 1/sum(w);
            avR = [avR Wmean]; avS = [avS Werr];
        end
    end
    % Fit we have more than 2 values get the chi^2
    num = length(avR);
    if num > 2
        savsta = [savsta sum(avR.^2./avS)/num];
        timsta = [timsta i*minav];
    end
end

% Fit the observed changes in dchi to an exponential function and estimate the asymptote
taus = [1 2 4 8 16 32 64 128 256 512 1024 2048 5096 10192];
ntau = length(taus);
nc = length(savsta);
for i = 1:ntau
    ef = (1-exp(-timsta/taus(i)));
    alpha(i) = mean(savsta./ef);
    res = savsta-(alpha(i)*ef);
    rmsa(i) = std(res);
end
[mnr, ir] = min(rmsa);
nrmsR = sqrt(alpha(ir))*dchi; tauf = taus(ir); 
fprintf(1,'NRMS Realistic %7.2f; Correlation time %8.2f days\n',nrmsR,taus(ir));

% % %
% ef = (1-exp(-timsta/taus(ir)));
% fit = ef*alpha(ir);
% 
%     figure;
%     plot(timsta,savsta,'rs',timsta,fit,'g^-');
%     axis tight, box on
%     xlabel('Averging time (days)');
%     ylabel('Chi^2 of Residuals');
%     set(gca,'Fontsize',15) 
%        xlabel('Averging time (days)')
%        ylabel('Chi^2 of Residuals')
end

        


