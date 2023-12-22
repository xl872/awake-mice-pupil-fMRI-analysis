clc;clear all;close all
V = BrikLoad ('ambmc_ori_100+orig'); 
LABEL = BrikLoad ('Allen_label_AMBMC_R+orig');
LABEL_L= BrikLoad ('Allen_label_AMBMC_L+orig'); 
LABEL_L(LABEL_L~=0)= LABEL_L(LABEL_L~=0)+865; 
LABEL=LABEL+LABEL_L;
V=reshape(mapminmax(V(:)',0,1),size(V));
az=45;el=45;
%%
load ADresultallnew.mat
Mcc_AD=Mcc;
cc_AD=cc;
load WTresultallnew.mat
Mcc_WT=Mcc;
cc_WT=cc;

% load ADresultrest.mat
% Mcc_AD=Mcc;
% cc_AD=cc;
% load WTresultrest.mat
% Mcc_WT=Mcc;
% cc_WT=cc;

% control
% CN=742;
% MO=12:16;
% ILA=190;
% KF=703;
% AUD=122:123;
% PERI=266:267;
% PH=626;
% AIP=221:224;
%  ROI={CN,MO,ILA,KF,AUD,PERI,PH,AIP,CN+865,MO+865,ILA+865,KF+865,AUD+865,PERI+865,PH+865,AIP+865};
VIS=160:166;%125:166;
RSP=232:250;
CA1=342:346;
PRNc=722;
CS=729:730;
MB=[725 734];%PRN,PRNr
LSV=444;
BST=477:492;
LGN=[516,552,553];
SC=[643,644];
ACA=175:176;
ROI={PRNc,PRNc+865,MB,MB+865,[SC,SC+865],VIS,VIS+865,[RSP,RSP+865],CA1,CA1+865,[ACA,ACA+865],LSV,LSV+865};

%all


% ROI=[ROI,ROI_c];
ROI_name={'VIS','RSP','PRNc','MB','CA1','LSV','SC'};
ROI_P=[];
for i=1:length(ROI)
    tROI=zeros(size(LABEL));
    for j=ROI{i}
        tROI(LABEL==j)=1;
    end
    ttROI=squeeze(sum(tROI,3));
    [a,b]=find(ttROI~=0);
    ROI_P(i,1)=mean(a);
    ROI_P(i,2)=mean(b);
    ttROI=squeeze(sum(tROI,[1 2]));
    c=find(ttROI);
    ROI_P(i,3)=mean(c);
end
% WT
pthr=0.01;
close all
figure;
subplot(1,3,1)
title 'WT'
VOX.IMG         = V;
VOX.ALPHA       = V.*0.035;
VOX.colormap    = gray;
% VOX_opts.ALPHA          = 0.3;
% VOX_opts.alpha_def      = 'standard';
VOX_opts.VoxelSize      = 1;

VOXview(VOX);
view(az,el);
camzoom(5)
axis off
hold on
linecolors = jet;
H=[];P=[];
ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
Mcc_n_WT=[];

for i=1:length(ROI)
    for j=1:length(ROI)
        % Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
        Mcc_n_WT(i,j)=round(mean(mean(Mcc_WT(ROI{i},ROI{j}),"omitnan"),"omitnan")*128)+128;
        %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
        [h(i,j),pWT(i,j)]=ttest(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2,"omitnan"),1,"omitnan")));
        if pWT(i,j)<pthr
            P=[ROI_P(i,:);ROI_P(j,:)];
            
            line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WT(i,j),:),'LineWidth',1.5)
            hold on
        end
    end
end

% AD

subplot(1,3,2)
title 'AD'
VOX.IMG         = V;
VOX.ALPHA       = V.*0.035;
VOX.colormap    = gray;
% VOX_opts.ALPHA          = 0.3;
% VOX_opts.alpha_def      = 'standard';
VOX_opts.VoxelSize      = 1;

VOXview(VOX);
view(az,el);
camzoom(5)
axis off
hold on
linecolors = jet;
H=[];P=[];
ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
Mcc_n_AD=[];

for i=1:length(ROI)
    for j=1:length(ROI)
        Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j}),"omitnan"),"omitnan")*128)+128;
        % Mcc_n_WT(i,j)=round(mean(mean(Mcc_WT(ROI{i},ROI{j})))*128)+128;
        %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
        [h(i,j),pAD(i,j)]=ttest(squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2,"omitnan"),1,"omitnan")));
        if pAD(i,j)<pthr
            P=[ROI_P(i,:);ROI_P(j,:)];
            
            line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_AD(i,j),:),'LineWidth',1.5)
            hold on
        end
    end
end


% WT-AD

subplot(1,3,3)
title 'WT - AD'
VOX.IMG         = V;
VOX.ALPHA       = V.*0.035;
VOX.colormap    = gray;
% VOX_opts.ALPHA          = 0.3;
% VOX_opts.alpha_def      = 'standard';
VOX_opts.VoxelSize      = 1;

VOXview(VOX);
view(az,el);
camzoom(5)
axis off
hold on
linecolors = jet;
H=[];P=[];
ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
Mcc_n_WTAD=[];

for i=1:length(ROI)
    for j=1:length(ROI)

        Mcc_n_WTAD(i,j)=round((mean(mean(Mcc_WT(ROI{i},ROI{j})))-mean(mean(Mcc_AD(ROI{i},ROI{j}))))*128)+128;
        [h(i,j),pWTAD(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));

        if pWTAD(i,j)<pthr
            P=[ROI_P(i,:);ROI_P(j,:)];
            
            line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WTAD(i,j),:),'LineWidth',1.5)
            hold on
        end
    end
end

% pWTAD=pWTAD(:);
% length(find(p<pthr));
%%
Mcc_n_WTAD=(Mcc_n_WTAD-128)/128;
Mcc_n_WT=(Mcc_n_WT-128)/128;
Mcc_n_AD=(Mcc_n_AD-128)/128;
%%
h=pWTAD;
h(h>pthr)=nan;
h(~isnan(h))=1;
figure;imagesc(Mcc_n_WTAD.*h);clim([-1 1]);colorbar;colormap('jet');%axis("xy")
%%
h=pWT;
h(h>pthr)=nan;
h(~isnan(h))=1;
figure;imagesc(Mcc_n_WT.*h);clim([-1 1]);colorbar;colormap('jet');%axis("xy")

h=pAD;
h(h>pthr)=nan;
h(~isnan(h))=1;
figure;imagesc(Mcc_n_AD.*h);clim([-1 1]);colorbar;colormap('jet');%axis("xy")


%%
% 
% load XL32.mat;
% 
% close all
% figure;
% subplot(1,3,1)
% title 'C1'
% VOX.IMG         = V;
% VOX.ALPHA       = V.*0.035;
% VOX.colormap    = gray;
% % VOX_opts.ALPHA          = 0.3;
% % VOX_opts.alpha_def      = 'standard';
% VOX_opts.VoxelSize      = 1;
% 
% VOXview(VOX);
% view(az,el);
% camzoom(5)
% axis off
% hold on
% linecolors = jet;
% H=[];P=[];
% ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
% ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
% Mcc_n_AD=[];Mcc_n_WT=[];
% XL1=reshape(mapminmax(XL(:,1)'),13,13);
% for i=1:length(ROI)
%     for j=1:length(ROI)
%         % Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
%         Mcc_n_WT(i,j)=round(XL1(i,j)*127)+128;
%         %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
% 
%             P=[ROI_P(i,:);ROI_P(j,:)];
% 
%             line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WT(i,j),:),'LineWidth',2)
%             hold on
% 
%     end
% end
% 
% subplot(1,3,2)
% title 'C2'
% VOX.IMG         = V;
% VOX.ALPHA       = V.*0.035;
% VOX.colormap    = gray;
% % VOX_opts.ALPHA          = 0.3;
% % VOX_opts.alpha_def      = 'standard';
% VOX_opts.VoxelSize      = 1;
% 
% VOXview(VOX);
% view(az,el);
% camzoom(5)
% axis off
% hold on
% linecolors = jet;
% H=[];P=[];
% ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
% ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
% Mcc_n_AD=[];Mcc_n_WT=[];
% XL2=reshape(mapminmax(XL(:,2)'),13,13);
% for i=1:length(ROI)
%     for j=1:length(ROI)
%         % Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
%         Mcc_n_WT(i,j)=round(XL2(i,j)*127)+128;
%         %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
% 
%             P=[ROI_P(i,:);ROI_P(j,:)];
% 
%             line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WT(i,j),:),'LineWidth',2)
%             hold on
% 
%     end
% end
% 
% subplot(1,3,3)
% title 'C3'
% VOX.IMG         = V;
% VOX.ALPHA       = V.*0.035;
% VOX.colormap    = gray;
% % VOX_opts.ALPHA          = 0.3;
% % VOX_opts.alpha_def      = 'standard';
% VOX_opts.VoxelSize      = 1;
% 
% VOXview(VOX);
% view(az,el);
% camzoom(5)
% axis off
% hold on
% linecolors = jet;
% H=[];P=[];
% ccn_AD=reshape(cell2mat(cc_AD),1730,1730,length(cc_AD));
% ccn_WT=reshape(cell2mat(cc_WT),1730,1730,length(cc_WT));
% Mcc_n_AD=[];Mcc_n_WT=[];
% XL3=reshape(mapminmax(XL(:,3)'),13,13);
% for i=1:length(ROI)
%     for j=1:length(ROI)
%         % Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
%         Mcc_n_WT(i,j)=round(XL3(i,j)*127)+128;
%         %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
% 
%             P=[ROI_P(i,:);ROI_P(j,:)];
% 
%             line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WT(i,j),:),'LineWidth',2)
%             hold on
% 
%     end
% end
