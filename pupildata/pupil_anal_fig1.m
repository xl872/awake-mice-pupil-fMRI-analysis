clc;clear all;close all;
load ADWT_fig1.mat

    %% model
    fps=30;
figure('Position',[100 100 500 150]);
subplot(1,2,1)
time=1:1200;time=time/30;


time_all=time;
ymean = AD.tm';
ystd= AD.ste;
plot(time_all,ymean,'-','linewidth',2,'Color',"b")
hold on
% h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
% set(h,'FaceColor',"b",'FaceAlpha',0.5,'EdgeColor','none');


ymean = WT.tm';
ystd= WT.ste;
plot(time_all,ymean,'-','linewidth',2,'Color',"r")
% 
% h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'r');
axis tight
ylim([0 1])
ylabel('Pupil size (a.u.)')
% set(h,'FaceColor',"r",'FaceAlpha',0.5,'EdgeColor','none');

hold off
legend('AD','WT',Location='southeast',Box='off',NumColumns=2)
xlabel('Time/s')




subplot(1,2,2)
tt=WT.tm;
keytime=WT.keytime;
keytime1=find(diff(tt)<0,1);
keytime2=find(tt(keytime1:end)<=1.25*min(tt(1:600)),1)+mean(keytime(:,:,1),[1 2])-1;
kt=find(zscore(tt(keytime2:end))>0.75*min(zscore(tt(keytime2:end))))+keytime2-1;
        keytime3=kt(find(kt>(keytime2+6*fps),1));
        kt=find(zscore(tt(keytime3:end))>0.75*max(zscore(tt(keytime3:end))))+keytime3-1;
        keytime4=kt(find(kt>(keytime3+0.2*fps),1));
        keytime4=1190;
ttt=1/fps:1/fps:1/fps*length(tt(keytime3:keytime4));%keytime(i,j,3)/fps:1/fps:40;%
 Ut=mean(tt(keytime4:keytime4+0.1*fps))%min(keytime(i,j,4)+10*fps,length(tt)))
        Lt=mean(tt(keytime3-0.1*fps:keytime3))

        ft = fittype(['b*(1-exp(-x/a))+',num2str(Lt)]);
        f=fit(ttt',smooth(tt(keytime3:keytime4),0.2),ft);
        p1=plot(f,'k:',ttt',tt(keytime3:keytime4),'r-'); 
        p1(1).LineWidth=2;p1(2).LineWidth=2;
        % legend('WT','fitted pupil dialation',Location='southeast',Box='off');
        axis tight
       axis([])
%         text(12,0.5,['{\itR}=',num2str(f.a)])

hold on
tt=AD.tm;
keytime=AD.keytime;
keytime1=find(diff(tt)<0,1);
keytime2=find(tt(keytime1:end)<=1.25*min(tt(1:600)),1)+mean(keytime(:,:,1),[1 2])-1;
kt=find(zscore(tt(keytime2:end))>0.75*min(zscore(tt(keytime2:end))))+keytime2-1;
        keytime3=kt(find(kt>(keytime2+6*fps),1));
        kt=find(zscore(tt(keytime3:end))>0.75*max(zscore(tt(keytime3:end))))+keytime3-1;
        % keytime4=kt(find(kt>(keytime3+0.2*fps),1));
        keytime4=1190;
ttt=1/fps:1/fps:1/fps*length(tt(keytime3:keytime4));%keytime(i,j,3)/fps:1/fps:40;%
 Ut=mean(tt(keytime4:keytime4+0.1*fps))%min(keytime(i,j,4)+10*fps,length(tt)))
        Lt=mean(tt(keytime3-0.1*fps:keytime3))

        ft = fittype(['b*(1-exp(-x/a))+',num2str(Lt)]);
        f=fit(ttt',smooth(tt(keytime3:keytime4),0.2),ft);
        p2=plot(f,'k:',ttt',tt(keytime3:keytime4),'b-');
        p2(1).LineWidth=2;p2(2).LineWidth=2;
      axis([])
%         text(9,0.5,['{\itR}=',num2str(f.a)])

     %% group analysis mean per mouse
     
figure('Position',[100 100 500 150]);
 subplot(1,3,2)
 per=75;
 rng('default')  % For reproducibility
    x1 = AD.cP;
    x2 = WT.cP;
    x1(x1<0)=nan;x2(x2<0)=nan;
    ii=0;x11=[];x22=[];
    for i = unique(ADindex)
        ii=ii+1;
        x11(ii) = trimmean(x1(:,ADindex==i),per,"all");
%         xt = harmmean(x1,1,"omitnan");
%         x11(ii) = harmmean(xt(ADindex==i),"all","omitnan");
    end
    ii=0;
    for i = unique(WTindex)
        ii=ii+1;
        x22(ii) = trimmean(x2(:,WTindex==i),per,"all");
%         xt = harmmean(x2,1,"omitnan");
%         x22(ii) = harmmean(xt(WTindex==i),"all","omitnan");
    end
    
    x = [x11'; x22'];
    g1 = repmat({'AD'},length(x11),1);
    g2 = repmat({'WT'},length(x22),1);
    g = [g1; g2];
   % boxplot(x,g)
   CategoricalScatterplot(x,g,'MarkerSize',8,'BoxLineWidth',1.5,...
   'MedianLineWidth',1.5,'WhiskerLineWidth',1.5)

   [h,p]= ttest2(x11,x22)
    title({'$P_d$'},'interpreter','latex')
    ylabel({'$P_d \rm{ (a.u.)}$'},'interpreter','latex')
ylim([0.1 0.8])

    subplot(1,3,3)
 rng('default')  % For reproducibility
    x1 = AD.R;
    x2 = WT.R;
    x1(x1>100)=nan;x2(x2>100)=nan;
    ii=0;x11=[];x22=[];
    for i = unique(ADindex)
        ii=ii+1;
        x11(ii) = trimmean(x1(:,ADindex==i),per,"all");
%         xt = harmmean(x1,1,"omitnan");
%         x11(ii) = harmmean(xt(ADindex==i),"all","omitnan");
    end
    ii=0;
    for i = unique(WTindex)
        ii=ii+1;
        x22(ii) = trimmean(x2(:,WTindex==i),per,"all");
%         xt = harmmean(x2,1,"omitnan");
%         x22(ii) = harmmean(xt(WTindex==i),"all","omitnan");
    end
    
    x = [x11'; x22'];
    g1 = repmat({'AD'},length(x11),1);
    g2 = repmat({'WT'},length(x22),1);

    g = [g1; g2];
    %boxplot(x,g)
    CategoricalScatterplot(x,g,'MarkerSize',8,'BoxLineWidth',1.5,...
   'MedianLineWidth',1.5,'WhiskerLineWidth',1.5)
    [h,p]= ttest2(x11,x22)
    %ylim([5 30])
       title({'$T$'},'interpreter','latex')
    ylabel({'$T \rm{ (a.u.)}$'},'interpreter','latex')

   subplot(1,3,1)
 rng('default')  % For reproducibility
    x1 = AD.dA;
    x2 = WT.dA;
    x1(x1<0)=nan;x2(x2<0)=nan;
    ii=0;x11=[];x22=[];
    for i = unique(ADindex)
        ii=ii+1;
        x11(ii) = trimmean(x1(:,ADindex==i),per,"all");
%         xt = harmmean(x1,1,"omitnan");
%         x11(ii) = harmmean(xt(ADindex==i),"all","omitnan");
    end
    ii=0;
    for i = unique(WTindex)
        ii=ii+1;
        x22(ii) = trimmean(x2(:,WTindex==i),per,"all");
%         xt = harmmean(x2,1,"omitnan");
%         x22(ii) = harmmean(xt(WTindex==i),"all","omitnan");
    end
    
    x = [x11'; x22'];
    g1 = repmat({'AD'},length(x11),1);
    g2 = repmat({'WT'},length(x22),1);
    g = [g1; g2];
    CategoricalScatterplot(x,g,'MarkerSize',8,'BoxLineWidth',1.5,...
   'MedianLineWidth',1.5,'WhiskerLineWidth',1.5)
    ylim([0.1 1])
    %boxplot(x,g)
    [h,p]= ttest2(x11,x22)
    title({'$P_c$'},'interpreter','latex')
    ylabel({'$P_c \rm{ (a.u.)}$'},'interpreter','latex')
    
