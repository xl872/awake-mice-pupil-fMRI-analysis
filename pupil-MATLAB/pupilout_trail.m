
clear all ;clc;close all;
%%
fps=30;
trail=[];
date='0810';EN=[22 23 24];mouse='2';%2 4
% trail=[trail;loadtrail(date,EN,fps)];
date='0810';EN=[32 33 34];mouse='4';%2 4
% trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[45 47 48]; mouse='7';%7 3 2
% trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[31 33 34]; mouse='3';%7 3 2
% trail=[trail;loadtrail(date,EN,fps)];
date='0811';EN=[16 17 18]; mouse='2';%7 3 2
% trail=[trail;loadtrail(date,EN,fps)];
date='0927';EN=[7 8 9]; mouse='2'; %2 3
% trail=[trail;loadtrail(date,EN,fps)];
date='0927';EN=[18 19 20]; mouse='3'; %2 3
% trail=[trail;loadtrail(date,EN,fps)];
date='1003';EN=[24 25]; mouse='3';
% trail=[trail;loadtrail(date,EN,fps)];
date='1025';EN=[12 13 14];mouse='6';
% trail=[trail;loadtrail(date,EN,fps)];
% date='1028';EN=[18 19 20];mouse='s4';
% % trail=[trail;loadtrail(date,EN,fps)];
name=date;
% name='WT';
output=mypupilout(date,date,EN,mouse);


%%
fps=30;
trail=[];
date='0729';EN=[13 14 15];mouse='9';
output=mypupilout(date,date,EN,mouse);
date='0731';EN=[10 11];mouse='7';
output=mypupilout(date,date,EN,mouse);

date='0810';EN=[22 23 24];mouse='2';%2 4
output=mypupilout(date,date,EN,mouse);
date='0810';EN=[32 33 34];mouse='4';%2 4
output=mypupilout(date,date,EN,mouse);
date='0811';EN=[45 47 48]; mouse='2';%7 3 2
output=mypupilout(date,date,EN,mouse);
date='0811';EN=[31 33 34]; mouse='3';%7 3 2
output=mypupilout(date,date,EN,mouse);
date='0811';EN=[16 17 18]; mouse='7';%7 3 2
output=mypupilout(date,date,EN,mouse);
date='0927';EN=[7 8 9]; mouse='2'; %2 3
output=mypupilout(date,date,EN,mouse);
date='0927';EN=[18 19 20]; mouse='3'; %2 3
output=mypupilout(date,date,EN,mouse);
date='1003';EN=[24 25]; mouse='3';
output=mypupilout(date,date,EN,mouse);
date='1025';EN=[12 13 14];mouse='6';
output=mypupilout(date,date,EN,mouse);
date='1028';EN=[18 19 20];mouse='s4';
output=mypupilout(date,date,EN,mouse);


%%
date='0111';EN=[19 20 21];mouse='7';
output=mypupilout(date,date,EN,mouse);
date='0125';EN=[10 19];mouse='12';
output=mypupilout(date,date,EN,mouse);
date='0125';EN=[30 31];mouse='13';
output=mypupilout(date,date,EN,mouse);
date='0309';EN=[15 16 17];mouse='7'; %
output=mypupilout(date,date,EN,mouse);
date='0309';EN=[26 27 28];mouse='8'; %
output=mypupilout(date,date,EN,mouse);
date='0310';EN=[33 34 35];mouse='8'; %
output=mypupilout(date,date,EN,mouse);
%
date='0330';EN=[26 27 28];mouse='9'; %
output=mypupilout(date,date,EN,mouse);
date='0330';EN=[36 37];mouse='4'; %
output=mypupilout(date,date,EN,mouse);
date='0331';EN=[15 16 17];mouse='9'; %
output=mypupilout(date,date,EN,mouse);
%%
fps=30;
trail=[];
date='0912';EN=[24 25 26];mouse='AD02';
output=mypupilout(date,date,EN,mouse);
date='0914';EN=[24 25]; mouse='AD03';% ad03 ad02
output=mypupilout(date,date,EN,mouse);
date='0914';EN=[43 45 46]; mouse='AD02';% ad03 ad02
output=mypupilout(date,date,EN,mouse);
date='0915';EN=[35 38 39 40]; mouse='AD03';
output=mypupilout(date,date,EN,mouse);

date='0926';EN=[10 11 12]; mouse='AD01';
output=mypupilout(date,date,EN,mouse);
date='1004';EN=[8 9 10]; mouse='AD04';
output=mypupilout(date,date,EN,mouse);
% date='1007';EN=[11 14]; mouse='AD02';
% output=mypupilout(date,date,EN,mouse);

date='1024';EN=[32 33 34];mouse='AD01';
output=mypupilout(date,date,EN,mouse);
date='1031';EN=[17 18 19];mouse='AD01'; %ad01
output=mypupilout(date,date,EN,mouse);
%date='1101';EN=[5 6 7];mouse='AD05'; %ad01
%output=mypupilout(date,date,EN,mouse);

%%
fps=30;
trail=[];
date='0223';EN=[6 7 8];mouse='ADN1'; %ad01
output=mypupilout(date,date,EN,mouse);
%
date='0224';EN=[15 16 17];mouse='AD02'; %
output=mypupilout(date,date,EN,mouse);
date='0224';EN=[46 47];mouse='ADN1'; %
output=mypupilout(date,date,EN,mouse);

date='0302';EN=[26 27];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);
date='0303';EN=[16 17 18];mouse='ADN5'; %
output=mypupilout(date,date,EN,mouse);
date='0303';EN=[27 28 29];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);
date='0309';EN=[38 39 40];mouse='ADN2'; %
output=mypupilout(date,date,EN,mouse);
date='0309';EN=[53 54 55];mouse='ADN3'; %
output=mypupilout(date,date,EN,mouse);
date='0310';EN=[10 11 12];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);

date='0330';EN=[9 10 11];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);
date='0331';EN=[5 6 7];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);

date='0407';EN=[13 14 15];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);

date='0421';EN=[34 35 36];mouse='ADN4'; %
output=mypupilout(date,date,EN,mouse);
%%
fps=30;
trail=[];
date='0912';EN=[24 25 26];mouse='AD02';
% trail=[trail;loadtrail(date,EN,fps)];
date='0914';EN=[24 25 43 45 46]; mouse='AD02';% ad03 ad02
% trail=[trail;loadtrail(date,EN,fps)];
date='0915';EN=[35 38 39 40]; mouse='AD03';
% trail=[trail;loadtrail(date,EN,fps)];
date='0926';EN=[10 11 12]; mouse='AD01';
% trail=[trail;loadtrail(date,EN,fps)];
date='1004';EN=[8 9 10]; mouse='AD04';
% trail=[trail;loadtrail(date,EN,fps)];
date='1007';EN=[11 12 14]; mouse='AD02';
% trail=[trail;loadtrail(date,EN,fps)];
date='1024';EN=[32 33 34];mouse='AD01';
% trail=[trail;loadtrail(date,EN,fps)];
date='1031';EN=[17 18 19];mouse='AD01'; %ad01
% trail=[trail;loadtrail(date,EN,fps)];
date='1101';EN=[5 6 7];mouse='AD05'; %ad01
% trail=[trail;loadtrail(date,EN,fps)];


% %
% name='AD';
name=date;
% mouse='ALL';
output=mypupilout(name,date,EN,mouse);
%%
function outputAM=mypupilout(name,date,EN,mouse)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% wt
year='2022';
fps=30;
i=1;
% date='0810';EN=[22 23 24];
% trail=[trail;loadtrail_all(date,EN,fps)];
for E = EN
trail=[];
trail=[trail;loadtrail_all(date,E,fps)];

% fps=30;
% i=1;trail=[];
% date='0927';EN=[18 19 20];
% trail=[trail;loadtrail_all(date,EN,fps)];
% date='1003';EN=[24 25];
% trail=[trail;loadtrail_all(date,EN,fps)];

% %% AD

% fps=30;
% i=1;trail=[];
% date='0912';EN=[24 25 26];
% trail=[trail;loadtrail_all(date,EN,fps)];
% date='0914';EN=[43 45 46];
% trail=[trail;loadtrail_all(date,EN,fps)];

% fps=30;
% i=1;trail=[];
% date='0926';EN=[10 11 12];
% trail=[trail;loadtrail_all(date,EN,fps)];
% date='1024';EN=[32 33 34];
% trail=[trail;loadtrail_all(date,EN,fps)];



for i=1
Ad1=diff(trail')';
Ad1=Ad1(:,10*fps:(20+40*10)*fps-1);
Ap=trail(:,10*fps+1:(20+40*10)*fps);


Ps=reshape(Ap,size(trail,1),60,205);Ps=squeeze(mean(Ps,2))';
Pd1=reshape(Ad1,size(trail,1),60,205);Pd1=squeeze(mean(Pd1,2))';

writematrix(Ps(:),['pupil.',mouse,'.',name,year,'E',num2str(E),'.dat'],'Delimiter',' ');  
% writematrix(Pd1(:),['pupild1.',name,'2022.dat'],'Delimiter',' ');  
end
%
pdc=[];
%    plot(Pd1(:,i));
pdc(:,1)=-Pd1([6,26,46,66,86,106,126,146,166,186]);

%
% WT only stim

fps=30;
i=1;trail=[];
% date='0810';EN=[22 23 24];
% trail=[trail;loadtrail(date,EN,fps)];


trail=[trail;loadtrail(date,E,fps)];


% fps=30;
% i=1;trail=[];
% date='0927';EN=[18 19 20];
% trail=[trail;loadtrail(date,EN,fps)];
% date='1003';EN=[24 25];
% trail=[trail;loadtrail(date,EN,fps)];

%AD only stim
% fps=30;
% i=1;trail=[];
% date='0912';EN=[24 25 26];
% trail=[trail;loadtrail(date,EN,fps)];
% date='0914';EN=[43 45 46];
% trail=[trail;loadtrail(date,EN,fps)];

% fps=30;
% i=1;trail=[];
% date='0926';EN=[10 11 12];
% trail=[trail;loadtrail(date,EN,fps)];
% date='1024';EN=[32 33 34];
% trail=[trail;loadtrail(date,EN,fps)];
%%
trail=mapminmax(trail,0,1);
M.trail=trail;

t1=reshape(trail,size(trail,1),40*fps,10);
tm=mean(mean(t1,3)',2);
t1=permute(t1,[3,1,2]);
t11=reshape(t1,size(t1,1)*size(t1,2),1200);
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
    end
end
%

T4=keytime(:,:,4)/fps;T2=keytime(:,:,2)/fps;
cP=(U-L);
M.c=c;M.cs=cs;M.lag=lag;M.lags=lags;M.keytime=keytime;M.L=L;M.U=U;M.cP=U-L;M.R=R;M.dA=dA;
%

%%
%%
for i=1:size(c,2)
   out_t= [mapminmax(1-c(:,i)',1,2);mapminmax(cP(:,i)',1,2);mapminmax(T4(:,i)'-T2(:,i)',1,2);mapminmax(abs(speed_c(:,i)'),1,2);mapminmax(abs(speed_d(:,i)'),1,2);mapminmax(U(:,i)',1,2);mapminmax(L(:,i)',1,2);mapminmax(pdc(:,i)',1,2)];
   output(i,:)=out_t(:);
   out_t= [mapminmax(dA(:,i)',1,2)]; 
   outputdA(i,:)=out_t(:);
   out_t= [mapminmax(cP(:,i)',1,2);mapminmax(R(:,i)',1,2)]; 
   outputAM(i,:)=out_t(:);
end
% 
fileID = fopen(['pupil_AM1dA.',mouse,'.',name,year,'E',num2str(E),'.txt'],'w');
fprintf(fileID,['6*%4.4f ' ...
    '46*%4.4f ' ...
    '86*%4.4f ' ...
    '126*%4.4f ' ...
    '166*%4.4f ' ...
    '206*%4.4f ' ...
    '246*%4.4f ' ...
    '286*%4.4f ' ...
    '326*%4.4f ' ...
    '366*%4.4f\n'], ...
    outputdA');
fclose(fileID);


fileID = fopen(['pupil_AM2fit.',mouse,'.',name,year,'E',num2str(E),'.txt'],'w');
fprintf(fileID,['6*%4.4f,%4.4f ' ...
    '46*%4.4f,%4.4f ' ...
    '86*%4.4f,%4.4f ' ...
    '126*%4.4f,%4.4f ' ...
    '166*%4.4f,%4.4f ' ...
    '206*%4.4f,%4.4f ' ...
    '246*%4.4f,%4.4f ' ...
    '286*%4.4f,%4.4f ' ...
    '326*%4.4f,%4.4f ' ...
    '366*%4.4f,%4.4f\n'], ...
    outputAM');
fclose(fileID);

end
end
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
%%
function trail=loadtrail(date,En,fps)
    i=1;trail=[];
    for en=En
        enum=['E',num2str(en)];
        load([date,enum,'_-21s.mat']);
        trail(i,:)=zscore(no_ev(Ps(20*fps+1:(20+40*10)*fps),0.005));
        i=i+1;
    end
end
function trail=loadtrail_205(date,En,fps)
    i=1;trail=[];
    for en=En
        enum=['E',num2str(en)];
        load([date,enum,'_-21s.mat']);
        trail(i,:)=zscore(no_ev(Ps(10*fps+1:(20+40*10)*fps),0.005));
        i=i+1;
    end
end

function trail=loadtrail_all(date,En,fps)
    i=1;trail=[];
    for en=En
        enum=['E',num2str(en)];
        load([date,enum,'_-21s.mat']);
        trail(i,:)=zscore(no_ev(Ps,0.005));
        i=i+1;
    end
end