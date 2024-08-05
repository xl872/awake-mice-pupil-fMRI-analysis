%% WT
clear all; close all;clc;
%%
fps=30;
trail=[];
date='0729';EN=[13 14 15];mouse='9';
trail=[trail;loadtrail(date,EN,fps)];
date='0731';EN=[10 11];mouse='7';
trail=[trail;loadtrail(date,EN,fps)];

date='0810';EN=[22 23 24];mouse='2';%2 4
trail=[trail;loadtrail(date,EN,fps)];
date='0810';EN=[32 33 34];mouse='4';%2 4
trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[45 47 48]; mouse='2';%7 3 2
trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[31 33 34]; mouse='3';%7 3 2
trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[16 17 18]; mouse='7';%7 3 2
trail=[trail;loadtrail(date,EN,fps)];
date='0927';EN=[7 8 9]; mouse='2'; %2 3
trail=[trail;loadtrail(date,EN,fps)];
date='0927';EN=[18 19 20]; mouse='3'; %2 3
trail=[trail;loadtrail(date,EN,fps)];
date='1003';EN=[24 25]; mouse='3';
trail=[trail;loadtrail(date,EN,fps)];
date='1025';EN=[12 13 14];mouse='6';
trail=[trail;loadtrail(date,EN,fps)];
date='1028';EN=[18 19 20];mouse='s4';
trail=[trail;loadtrail(date,EN,fps)];
date='0111';EN=[19 20 21];mouse='7';
trail=[trail;loadtrail(date,EN,fps)];
date='0125';EN=[10 19];mouse='s12';%10 15 19
trail=[trail;loadtrail(date,EN,fps)];
date='0125';EN=[30 31];mouse='s13';
trail=[trail;loadtrail(date,EN,fps)];
date='0309';EN=[15 16 17];mouse='7'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0309';EN=[26 27 28];mouse='8'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0310';EN=[33 34 35];mouse='8'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0330';EN=[26 27 28];mouse='9'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0330';EN=[36 37];mouse='4'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0331';EN=[15 16 17];mouse='9'; %
trail=[trail;loadtrail(date,EN,fps)];
name=date;
name='WT';
trail=mapminmax(trail,0,1);
M.trail=trail;
trail=trail(:,[10*fps+1:410*fps]);

t1=reshape(trail,size(trail,1),40*fps,10);
tm=mean(mean(t1,3)',2);
t1=permute(t1,[3,1,2]);
t11=reshape(t1,size(t1,1)*size(t1,2),1200);
M.tm=tm;
M.ste=std(t11)/square(size(t11,1));
% corr in pupil(for stats)
s=[zeros(8*30,1);ones(32*30,1)];
c=[];lag=[];R=[];dA=[];
cs=[];lags=[];keytime=[];

% figure;
% time=1:1200;time=time/30;


for i=1:size(t1,1)
    for j=1:size(t1,2)
        tt=squeeze(t1(i,j,:));
        [tc,tlag] = xcorr(tm',tt','normalized');
        [tcs,tlags] = xcorr(s',tt','normalized');
        c(i,j)=max(tc);
        cs(i,j)=max(tcs);
        lag(i,j)=tlag(find(tc==max(tc),1))/30;
        lags(i,j)=tlags(find(tcs==max(tcs),1))/30;
        tt=smooth(zscore(squeeze(t1(i,j,:))));
        keytime(i,j,1)=find(diff(smooth(tt))<0,1);
        keytime(i,j,2)=find(tt(keytime(i,j,1):end)<=0.75*min(tt(1:600)),1)+keytime(i,j,1)-1;
        kt=find(zscore(tt(keytime(i,j,2):end))>0.75*min(zscore(tt(keytime(i,j,2):end))))+keytime(i,j,2)-1;
        keytime(i,j,3)=kt(find(kt>(keytime(i,j,2)+2*fps),1));
        kt=find(zscore(tt(keytime(i,j,3):end))>0.75*max(zscore(tt(keytime(i,j,3):end))))+keytime(i,j,3)-1;
        keytime(i,j,4)=kt(find(kt>(keytime(i,j,3)+0.2*fps),1));
        time_c(i,j)=keytime(i,j,2)-keytime(i,j,1);
        time_s(i,j)=keytime(i,j,3)-keytime(i,j,2);
        time_d(i,j)=keytime(i,j,4)-keytime(i,j,3);
        speed=diff(tt);
        speed_d(i,j)=mean(speed(keytime(i,j,3):keytime(i,j,4)));
        speed_c(i,j)=mean(speed(keytime(i,j,1):keytime(i,j,2)));
        speed_dm(i,j)=max(speed(keytime(i,j,3):keytime(i,j,4)));
        speed_cm(i,j)=min(speed(keytime(i,j,1):keytime(i,j,2)));
        tt=smooth(squeeze(t1(i,j,:)));
        U(i,j)=mean(tt(keytime(i,j,4):keytime(i,j,4)+3));%min(keytime(i,j,4)+10*fps,length(tt)))
        L(i,j)=mean(tt(keytime(i,j,3)-3:keytime(i,j,3)));
        ttt=1/fps:1/fps:1/fps*length(tt(keytime(i,j,3):keytime(i,j,4)));%keytime(i,j,3)/fps:1/fps:40;%
        ft = fittype(['b*(1-exp(-x/a))+',num2str(L(i,j))]);
        f=fit(ttt',smooth(tt(keytime(i,j,3):keytime(i,j,4)),0.2),ft,'Start', [11,U(i,j)-L(i,j)],'Lower',[0 0],'Robust', 'LAR');
%         figure;plot(f,ttt',tt(keytime(i,j,3):min(keytime(i,j,4)+4*fps,length(tt))));
        R(i,j)=f.a;
        dP=f.b;
        dA(i,j)=(mean(tt(1:keytime(i,j,1)))-mean(L(i,j)))/mean(tt(1:keytime(i,j,1)));
        
%         plot(time,squeeze(t1(i,j,:)),LineWidth=0.5,Color=[1 0.75 0.75 0.5]); hold on
        
%         sortt=sort(tt);
%         L(i,j)=mean(sortt(end-119:end));
%         U(i,j)=mean(sortt(1:120));
    end
end
%

T4=keytime(:,:,4)/fps;T2=keytime(:,:,2)/fps;
cP=(U-L);
M.c=c;M.cs=cs;M.lag=lag;M.lags=lags;M.keytime=keytime;M.L=L;M.U=U;M.cP=U-L;M.R=R;M.dA=dA;M.dP=dP;M.sdm=speed_dm;M.scm=speed_cm;
%
%
WT=M;
%%
fps=30;
trail=[];
date='0912';EN=[24 25 26];mouse='AD02';
trail=[trail;loadtrail(date,EN,fps)];
date='0914';EN=[24 25]; mouse='AD03';% ad03 ad02
trail=[trail;loadtrail(date,EN,fps)];
date='0914';EN=[43 45 46]; mouse='AD02';% ad03 ad02
trail=[trail;loadtrail(date,EN,fps)];
date='0915';EN=[35 38 39 40]; mouse='AD03';
trail=[trail;loadtrail(date,EN,fps)];
date='0926';EN=[10 11 12]; mouse='AD01';
trail=[trail;loadtrail(date,EN,fps)];
date='1004';EN=[8 9 10]; mouse='AD04';
trail=[trail;loadtrail(date,EN,fps)];
% date='1007';EN=[11 12 14]; mouse='AD02';
% trail=[trail;loadtrail(date,EN,fps)];
date='1024';EN=[32 33 ];mouse='AD01';%34
trail=[trail;loadtrail(date,EN,fps)];
date='1031';EN=[17 18 19];mouse='AD01'; %ad01
trail=[trail;loadtrail(date,EN,fps)];

date='0223';EN=[6 7 8];mouse='ADN1'; %ad01
trail=[trail;loadtrail(date,EN,fps)];
%
date='0224';EN=[15 16 17];mouse='AD02'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0224';EN=[46 47];mouse='ADN1'; %
trail=[trail;loadtrail(date,EN,fps)];

date='0302';EN=[26 27];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0303';EN=[16 17 18];mouse='ADN5'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0303';EN=[27 28 29];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];

date='0309';EN=[38 39 40];mouse='ADN2'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0309';EN=[53 54 55];mouse='ADN3'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0310';EN=[10 11 12];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0330';EN=[9 10 11];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0331';EN=[5 6 7];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];
date='0407';EN=[13 14 15];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];

date='0421';EN=[34 35 36];mouse='ADN4'; %
trail=[trail;loadtrail(date,EN,fps)];
name='AD';
% mouse='ALL';
trail=mapminmax(trail,0,1);
M.trail=trail;
trail=trail(:,[10*fps+1:410*fps]);

t1=reshape(trail,size(trail,1),40*fps,10);
tm=mean(mean(t1,3)',2);
t1=permute(t1,[3,1,2]);
t11=reshape(t1,size(t1,1)*size(t1,2),1200);
M.tm=tm;
M.ste=std(t11)/square(size(t11,1));

% corr in pupil(for stats)
s=[zeros(8*30,1);ones(32*30,1)];
c=[];lag=[];R=[];dA=[];
cs=[];lags=[];keytime=[];
for i=1:size(t1,1)
    for j=1:size(t1,2)
        tt=squeeze(t1(i,j,:));
        [tc,tlag] = xcorr(tm',tt','normalized');
        [tcs,tlags] = xcorr(s',tt','normalized');
        c(i,j)=max(tc);
        cs(i,j)=max(tcs);
        lag(i,j)=tlag(find(tc==max(tc),1))/30;
        lags(i,j)=tlags(find(tcs==max(tcs),1))/30;
        tt=smooth(zscore(squeeze(t1(i,j,:))));
        keytime(i,j,1)=find(diff(smooth(tt))<0,1);
        if find(tt(keytime(i,j,1):end)<=0.75*min(tt(1:600)),1)
        keytime(i,j,2)=find(tt(keytime(i,j,1):end)<=0.75*min(tt(1:600)),1)+keytime(i,j,1)-1;
        else
            keytime(i,j,2)=keytime(i,j,1)+8*30;
        end
        kt=find(zscore(tt(keytime(i,j,2):end))>0.75*min(zscore(tt(keytime(i,j,2):end))))+keytime(i,j,2)-1;
        keytime(i,j,3)=kt(find(kt>(keytime(i,j,2)+3*fps),1));
        kt=find(zscore(tt(keytime(i,j,3):end))>0.75*max(zscore(tt(keytime(i,j,3):end))))+keytime(i,j,3)-1;
        keytime(i,j,4)=kt(find(kt>(keytime(i,j,3)+0.2*fps),1));
        time_c(i,j)=keytime(i,j,2)-keytime(i,j,1);
        time_s(i,j)=keytime(i,j,3)-keytime(i,j,2);
        time_d(i,j)=keytime(i,j,4)-keytime(i,j,3);
        speed=diff(tt);
        speed_d(i,j)=mean(speed(keytime(i,j,3):keytime(i,j,4)));
        speed_c(i,j)=mean(speed(keytime(i,j,1):keytime(i,j,2)));
        speed_dm(i,j)=max(speed(keytime(i,j,3):keytime(i,j,4)));
        speed_cm(i,j)=min(speed(keytime(i,j,1):keytime(i,j,2)));
        tt=smooth(squeeze(t1(i,j,:)));
        U(i,j)=mean(tt(keytime(i,j,4):max((keytime(i,j,4)+3),length(tt))));%min(keytime(i,j,4)+10*fps,length(tt)))
        L(i,j)=mean(tt(keytime(i,j,3)-3:keytime(i,j,3)));
        ttt=1/fps:1/fps:1/fps*length(tt(keytime(i,j,3):keytime(i,j,4)));%keytime(i,j,3)/fps:1/fps:40;%
        ft = fittype(['b*(1-exp(-x/a))+',num2str(L(i,j))]);
        f=fit(ttt',smooth(tt(keytime(i,j,3):keytime(i,j,4)),0.2),ft,'Start', [11,U(i,j)-L(i,j)],'Lower',[0 0],'Robust', 'LAR');
%         figure;plot(f,ttt',tt(keytime(i,j,3):min(keytime(i,j,4)+4*fps,length(tt))));
        R(i,j)=f.a;
        dP=f.b;
        dA(i,j)=(mean(tt(1:keytime(i,j,1)))-mean(L(i,j)))/mean(tt(1:keytime(i,j,1)));
%         sortt=sort(tt);
%         L(i,j)=mean(sortt(end-119:end));
%         U(i,j)=mean(sortt(1:120));
%         plot(time,squeeze(t1(i,j,:)),LineWidth=0.5,Color=[0.7 0.8 0.9 0.5]); hold on
        
    end
end
%

T4=keytime(:,:,4)/fps;T2=keytime(:,:,2)/fps;
cP=(U-L);
M.c=c;M.cs=cs;M.lag=lag;M.lags=lags;M.keytime=keytime;M.L=L;M.U=U;M.cP=U-L;M.R=R;M.dA=dA;M.dP=dP;M.sdm=speed_dm;M.scm=speed_cm;
%
%
AD=M;

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
    
 



%%
function [x1] = no_ev(x,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x1=x;
x=sort(x);
n=round(p*length(x));
minx=x(n+1);
maxx=x(length(x)-n);
x1(x1>maxx)=maxx;x1(x1<minx)=minx;
end

function trail=loadtrail(date,En,fps)
    i=1;trail=[];
    for en=En
        enum=['E',num2str(en)];
        load([date,enum,'_-21s.mat']);
        trail(i,:)=zscore(no_ev(Ps(9*fps:(19+40*10)*fps-1),0.001));
        i=i+1;
    end
end

