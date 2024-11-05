% Plot joint posterior pdfs

l=length(MCMC(1,:));
step=1;
pl=ceil(l/step);
Cmax = [];
for i = 1:step:l
    for j = 1:step:l
        if i==j
            subplot(pl,pl,pl*(ceil(i/step)-1)+ceil(j/step))
            h1=histogram(MCMC(:,i),'normalization','probability','Linestyle','none','Linewidth',0.1,'Displaystyle','bar','Edgecolor','b'); axis tight, hold on
            box on, grid on
            dummy=h1.BinEdges(h1.Values==max(h1.Values));
            Cmax(i,1)=dummy(1);
            Cmax(i,[2 3]) = prctile(MCMC(:,i),[2.5 97.5]);
            plot([Cmax(i) Cmax(i)],get(gca,'ylim'),'r','lineWidth',1.5);
            if i==l
                xlabel(names(i))
                set(gca,'YTickLabel','','FontSize',15)
            else
                set(gca,'XTickLabel','','YTickLabel','')
            end            
            ylim([0 max(h1.Values)])
            xlim(h1.BinLimits)
        elseif i>j
            subplot(pl,pl,pl*(ceil(i/step)-1)+ceil(j/step))
            h=histogram2(MCMC(:,j),MCMC(:,i),'Normalization','pdf','FaceColor','flat','LineStyle','none','DisplayStyle','tile','ShowEmptyBins','off'); hold on
            if i==l
                xlabel(names(j))
                set(gca,'YTickLabel','','FontSize',15)
            else
                set(gca,'XTickLabel','','YTickLabel','')
            end
            
            if j==1
                ylabel(names(i))
                set(gca,'FontSize',15,'YTickLabelMode','auto')
            end
            %colormap hot
            colormap(flipud(ttscm('davos')))
            axis tight
            grid off
        else
        end
    end
end