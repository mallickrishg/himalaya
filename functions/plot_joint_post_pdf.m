function [msol] = plot_joint_post_pdf(MCMC,names)
% Function to plot 1-d and 2-d marginal posterior pdfs, and returns MAP 
% INPUTS
% MCMC : is a matrix of size Nsamples X Nparams
% names : cell array containing names of all Nparams
% OUTPUT - lower triangle matrix of figures has 2-D marginals, diagonal
% contains 1-D marginal odfs
% msol - MAP and 95%CI of model parameters (Nparams x 3) - column1 - MAP,
% column2-3 95%CI
% Rishav Mallick, EOS 2018

if isempty(names)
    for i = 1:length(MCMC(1,:))
        names{i} = num2str(i);
    end
end

l=length(MCMC(1,:));
step=1;
pl=ceil(l/step);
Cmax = [];

[ha,~] = tight_subplot(pl,pl,0.05,[0.1 0.01],0.1);
for i = 1:step:l
    for j = 1:step:l
        if i==j
            %subplot(pl,pl,pl*(ceil(i/step)-1)+ceil(j/step))
            axes(ha(pl*(ceil(i/step)-1)+ceil(j/step)))
            
            h1=histogram(MCMC(:,i),'normalization','probability','Linestyle','none','Linewidth',0.1,'Displaystyle','bar','Edgecolor','b'); axis tight, hold on
            box on, %grid on
            dummy=h1.BinEdges(h1.Values==max(h1.Values)) + h1.BinWidth/2;
            Cmax(i,1)=dummy(1);
            Cmax(i,[2 3]) = prctile(MCMC(:,i),[2.5 97.5]);
            plot([Cmax(i) Cmax(i)],get(gca,'ylim'),'r','lineWidth',1.5);
            if i==l
                xlabel(names(i))
                set(gca,'YTickLabel','','FontSize',12)
            else
                set(gca,'XTickLabel','','YTickLabel','')
            end            
            ylim([0 max(h1.Values)])
            xlim(h1.BinLimits)
        elseif i>j
            %subplot(pl,pl,pl*(ceil(i/step)-1)+ceil(j/step))
            axes(ha(pl*(ceil(i/step)-1)+ceil(j/step)))
            
            h=histogram2(MCMC(:,j),MCMC(:,i),'Normalization','pdf','FaceColor','flat','LineStyle','none','DisplayStyle','tile','ShowEmptyBins','off'); hold on
            if i==l
                xlabel(names(j))
                set(gca,'YTickLabel','','FontSize',12)
            else
                set(gca,'XTickLabel','','YTickLabel','')
            end
            
            if j==1
                ylabel(names(i))
                set(gca,'FontSize',12,'YTickLabelMode','auto')
            end
            %colormap hot
            %colormap(flipud(ttscm('davos')))
            axis tight
            grid off
        elseif i<j
            ax = (ha(pl*(ceil(i/step)-1)+ceil(j/step)));
            ax.Visible = 'off';
        end
    end
end

msol = Cmax;
end