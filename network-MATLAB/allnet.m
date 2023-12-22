clc;clear all;close all
V = BrikLoad ('ambmc_ori_100+orig'); 
LABEL = BrikLoad ('Allen_label_AMBMC_R+orig');
LABEL_L= BrikLoad ('Allen_label_AMBMC_L+orig'); 
LABEL_L(LABEL_L~=0)= LABEL_L(LABEL_L~=0)+865; 
LABEL=LABEL+LABEL_L;
V=reshape(mapminmax(V(:)',0,1),size(V));
az=90;el=90;

T = readtable('labels.txt');
%
labelnew={};
par=digitsPattern;
for i=3:866
    labelt=table2array(T(i,2));
    ni=strfind(labelt{:},par);
    if isempty(ni)
        labelnew(i)={labelt{:}};
    else
        labelnew(i)={labelt{:}(1:ni(1)-1)};
    end
end

%
% labelnew={};
% par=digitsPattern;
% for i=5:866
%     labelt=table2array(T(i,2));
%     ni1=strfind(labelt{:},par);
%     labelt=labelt{:};
%     ni2=find(labelt<='z'&labelt>='a');
%     ni=sort([ni1,ni2]);
%     if isempty(ni)
%         labelnew(i)={labelt};
%     else
%         labelnew(i)={labelt(1:ni(1)-1)};
%     end
% end
%%
ROI={};
ROInow='';
ROIn=0;
for i=9:866
    if ~strcmp(labelnew{i},ROInow)
        ROIn=ROIn+1;
        ROI(ROIn,1)=labelnew(i);
        ROI(ROIn,2)={i-1};
        ROInow=labelnew{i};
    else
        ROI(ROIn,2)={cat(2,cell2mat(ROI(ROIn,2)),i-1)};
    end
end

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


ROI_my=ROI;
ROI=[ROI(:,2)',cellfun(@(x) x + 865, ROI(:,2)', 'UniformOutput', false)];

ROI_name=ROI_my(:,1)';
%%
athr=11;
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
%%
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
Mcc_n_AD=[];Mcc_n_WT=[];
ROI_all=[ROI_my;ROI_my];
ROI_S={};count=0;
for i=1:length(ROI)
    for j=1:length(ROI)
        % Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
        Mcc_n_WT(i,j)=round(mean(mean(Mcc_WT(ROI{i},ROI{j}),"omitnan"),"omitnan")*128)+128;
        %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
        [h(i,j),p(i,j)]=ttest(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)));
    end
end
asort=sort(Mcc_n_WT(:));
for i=1:length(ROI)
    for j=1:length(ROI)
      
        if p(i,j)<pthr & (Mcc_n_WT(i,j)>asort(end-athr) | Mcc_n_WT(i,j)<asort(athr+1))
            P=[ROI_P(i,:);ROI_P(j,:)];
            count=count+1;
            ROI_S(count,1)=ROI_all(i,1);
            ROI_S(count,2)=ROI_all(j,1);
            line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WT(i,j),:),'LineWidth',2)
            hold on
        end
    end
end

ROI_SWT=ROI_S;
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
Mcc_n_AD=[];Mcc_n_WT=[];
ROI_S={};count=0;
for i=1:length(ROI)
    for j=1:length(ROI)
        Mcc_n_AD(i,j)=round(mean(mean(Mcc_AD(ROI{i},ROI{j})))*128)+128;
        % Mcc_n_WT(i,j)=round(mean(mean(Mcc_WT(ROI{i},ROI{j})))*128)+128;
        %[H(i,j),P(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
        [h(i,j),p(i,j)]=ttest(squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));
    end
end
asort=sort(Mcc_n_AD(:));
for i=1:length(ROI)
    for j=1:length(ROI)
        if p(i,j)<pthr & (Mcc_n_AD(i,j)>asort(end-athr) | Mcc_n_AD(i,j)<asort(athr+1))
            P=[ROI_P(i,:);ROI_P(j,:)];
            count=count+1;
            ROI_S(count,1)=ROI_all(i,1);
            ROI_S(count,2)=ROI_all(j,1);
            line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_AD(i,j),:),'LineWidth',2)
            hold on
        end
    end
end
ROI_SAD=ROI_S;

 %% WT-AD abs
 close all
athr=200*2;
pthr=0.01;
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
Mcc_n_AD=[];Mcc_n_WT=[];Mcc_n_WTAD=[];
ROI_S={};count=0;
for i=1:length(ROI)
    for j=1:length(ROI)

        Mcc_n_WTAD(i,j)=abs((mean(mean(Mcc_WT(ROI{i},ROI{j})))-mean(mean(Mcc_AD(ROI{i},ROI{j})))));
        [h(i,j),p(i,j)]=ttest2(squeeze(mean(mean(ccn_WT(ROI{i},ROI{j},:),2),1)),squeeze(mean(mean(ccn_AD(ROI{i},ROI{j},:),2),1)));

    end
end
%
A=zeros(size(p));
A(p<pthr)=1;
asort=sort(Mcc_n_WTAD(:).*A(:));
for i=1:length(ROI)
    for j=1:length(ROI)
        if p(i,j)<pthr & (Mcc_n_WTAD(i,j)>asort(end-athr))% | Mcc_n_WTAD(i,j)<asort(athr+1))
            P=[ROI_P(i,:);ROI_P(j,:)];
            count=count+1;
            ROI_S(count,1)=ROI_all(i,1);
            ROI_S(count,2)=ROI_all(j,1);
            % line(P(:,2),P(:,1),P(:,3),'Color', linecolors(Mcc_n_WTAD(i,j),:),'LineWidth',2)
            hold on
        end
    end
end
ROI_SWTAD=ROI_S;

%

ROI_out={};
ROInow='';
ROIn=0;
ROI_SSWTAD=sort(ROI_SWTAD(:,1));
for i=1:length(ROI_SWTAD)
    if ~strcmp(ROI_SSWTAD{i},ROInow)
        ROIn=ROIn+1;
        ROI_out(ROIn,1)=ROI_SSWTAD(i);
        ROI_out(ROIn,2)={1};
        ROInow=ROI_SSWTAD{i};
    else
        ROI_out(ROIn,2)={ROI_out{ROIn,2}+1};
    end
end
%
figure;pie(cell2mat(ROI_out(:,2)))
 legend(ROI_out(:,1),'Box','off');
%  fontsize(20,"points")
 %
 figure;pie(categorical(ROI_SSWTAD),{'PRNc','PRNr'})
 % fontsize(20,"points")