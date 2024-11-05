function [model,Rco,phdel] = calc_model_GRACE_TRMM(n,isort,ind,maxrowplot,maxcolplot,tvec_gr,dm,ts_TRMM)
% calc_model_GRACE_TRMM can be used to calculate a simple linear least-squares model to
% GRACE & TRMM data and plot them based on latitude and longitude
% INPUTS - 
% n - total number of stations/grid points
% isort - sorted array of lat,lon. Here, I use latitude (isort(:,1)) sort to plot
% ind - sorted indices
% maxrow,maxcol - number of rows/columns for subplot
% tvec_gr - cell array of time-vectors (julian-day) for GRACE data
% dm - cell array of mass change (GT) for GRACE data
% ts_TRMM - array of structures that contains tvec
% (julian-day),surface-precipitation data (mm), lonGR,latGR - location at
% which data was sampled
% OUTPUTS -
% model - structure containing tvec,predicted - rfall,dm and model
% parameters Mdm [ offset, slope, A1, A2], Mrfall [mean, A1 A2 B1 B2]
% Rishav Mallick, EOS, 2018

% define anonymous functions to create design matrix
G = @(t) [ones(length(t),1) t sin(2*pi*t/365) cos(2*pi*t/365)];
G2 = @(t) [ones(length(t),1) sin(2*pi*t/365) cos(2*pi*t/365) sin(4*pi*t/365) cos(4*pi*t/365)];
% GigaTons to m
GT2m = .080623;

% declare variables to calculate some simple statistics
% phase delay between masschange and rfall
phdel = [];
maxdel = 6; % 6 months - maximum delay
Rco = []; % corr-coeff

close all
figure(1),clf
set(gcf,'position',get(0,'Screensize'))

for i = 1:n
    % extract minmax times of overlapping datasets
    tmin = max([min(tvec_gr{ind(i)}) min(ts_TRMM(ind(i)).tvec)]);
    tmax = min([max(tvec_gr{ind(i)}) max(ts_TRMM(ind(i)).tvec)]);
    tminindex = find(abs(tvec_gr{ind(i)}-tmin)<=1,1,'first');
    tmaxindex = find(abs(tvec_gr{ind(i)}-tmax)<=1,1,'first');
    
    tplot = tvec_gr{ind(i)}; tplot = tplot(tminindex:tmaxindex);
    
    % extract appropriate data and interpolate to get comparable data
    dmplot = dm{ind(i)}.*GT2m; dmplot = dmplot(tminindex:tmaxindex);
    rfall = interp1(ts_TRMM(ind(i)).tvec,ts_TRMM(ind(i)).data./1e3,tplot);
    
    Mdm = (G(tplot)\dmplot);
    preddm = G(tplot)*Mdm;
    
    Mrfall = (G2(tplot)\rfall);
    predrfall = G2(tplot)*Mrfall;
    
    model(i).rfall = predrfall;
    model(i).dm = preddm;
    model(i).tvec = tplot;
    model(i).Mdm = Mdm;
    model(i).Mrfall = Mrfall;
    
    %calulate phase delay in signals, and corr-coefficient
    phdel(i,1) = finddelay(dmplot,rfall,maxdel);
    Rco(i,1) = corr(dmplot,rfall);
    
    if i==1
        subplot(maxrowplot,maxcolplot,i)
        
        plot(tplot,dmplot,'k','LineWidth',.5), hold on, axis tight
        plot(tplot,preddm,'r','LineWidth',.1)
        ylabel('\Delta m (m)'), ylim([-.4 .4])
        
        yyaxis right
        plot(tplot,rfall,'Color',rgb('royalblue'),'LineWidth',.5), hold on
        plot(tplot,predrfall,'-','Color',rgb('forestgreen'),'LineWidth',.1)
        datetick('x','yy'), axis tight
        ylim([0 1.5]), ylabel('Precip (m)')
        
        plotnum = 2;
    
    elseif isort(i,1)==isort(i-1,1)
        subplot(maxrowplot,maxcolplot,plotnum)
        
        plot(tplot,dmplot,'k','LineWidth',.5), hold on, axis tight
        plot(tplot,preddm,'r','LineWidth',.1)
        ylim([-.4 .4]),ylabel('')
        set(gca,'YTickLabel','')
        
        yyaxis right
        plot(tplot,rfall,'Color',rgb('royalblue'),'LineWidth',.5), hold on
        plot(tplot,predrfall,'-','Color',rgb('forestgreen'),'LineWidth',.1)
        datetick('x','yy'), axis tight, box on
        ylim([0 1.5])
        set(gca,'YTickLabel','')
        
        plotnum = plotnum+1;
    
    else
        k = ceil(plotnum/maxcolplot);
        if mod(plotnum,maxcolplot)~=1   
            plotnum = k*maxcolplot +1;
        end
        subplot(maxrowplot,maxcolplot,plotnum)
        
        plot(tplot,dmplot,'k','LineWidth',.5), hold on
        plot(tplot,preddm,'r','LineWidth',.1)
        datetick('x','yy'), axis tight
        ylabel('\Deltam (m)'), ylim([-.4 .4])
        
        yyaxis right
        plot(tplot,rfall,'Color',rgb('royalblue'),'LineWidth',.5), hold on
        plot(tplot,predrfall,'-','Color',rgb('forestgreen'),'LineWidth',.1)
        datetick('x','yy'), axis tight
        ylim([0 1.5]), ylabel('Precip (m)')
        
        plotnum = plotnum + 1;
    end
    title(['(' num2str(ts_TRMM(ind(i)).lonGR) ',' num2str(ts_TRMM(ind(i)).latGR) ')'],'FontWeight','normal')
    set(gca,'FontSize',8,'XTickLabelRotation',45,'Color','none')
end

end