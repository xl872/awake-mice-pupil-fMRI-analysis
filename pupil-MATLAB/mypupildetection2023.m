% Before detection, have a look at the pupil image and find the pupil/reference area manually.
close all; clear all; clc;
 %load the movie  
date='0912';
enum='E24'; 
 vd=['D:/22/U19_09022022/scale/videoAD/',date,'2022dummmy=40_',enum,'.avi'];
rotate=45;
v = VideoReader(vd);
f=[];frame=[];
for i=1:217*30
    f = uint8(imrotate(readFrame(v),rotate,'bilinear','crop'));
    
end
frame=f;
rp=squeeze(mean(mean(frame(:,:,1,:))));
figure;image(f)

%% detection
close all; clear all; clc;

for en=[24]
%load the movie & have a look 
enum=['E',num2str(en)]
    rotate=45;
date='0912';
 vd=['D:/22/U19_09022022/scale/videoAD/',date,'2022dummmy=40_',enum,'.avi'];

v = VideoReader(vd);
f=[];frame=[];
for i=1:8*30 %chose the frame. the 8 second and 30 is the sample rate of the camera
    f = uint8(imrotate(readFrame(v),rotate,'bilinear','crop'));
    frame(:,:,:,i)=f;
end
rp=squeeze(mean(mean(frame(:,:,1,:))));
figure;image(f)

%% set the area

% eq='f(225:284,316:417,:)'; %pupil area begin
% eqb='f(228:303,63:303,'; %reference area begin
% eq='f(227:292,307:408,:)';%end
% eqb='f(279:303,138:427,';

eq='f(216:299,319:419,:)';
eqb='f(228:303,63:303,'; 
% eqb='f(:,:,'; %reference area
figure;eval(['image(',eq,')']);

%%
ctst=[0.001 0.12];%ctst=[0.05 0.2]; contrast threshold setting of pupil area
ctstr=[0.01 0.3]; %contrast threshold setting of reference area


v = VideoReader(vd);
i=1;frame=[];rp=[];r=[];ap=[];a=[];gp=[];g=[];b=[];bp=[];
while hasFrame(v) & i<=60*30*8 %8 mins  
    f = uint8(imrotate(readFrame(v),rotate,'bilinear','crop'));
    eval(['frame(:,:,:,i)=',eq,';']);    
    eval(['tf=',eq,';']);
    ga=rgb2gray(tf);
%seperate the channels 
%pupil
    ap(i)=mean(mean(imadjust(ga,ctst*2)));
    rp(i)=mean(mean(imadjust(tf(:,:,1),ctst*3)));
    gp(i)=mean(mean(imadjust(tf(:,:,2),ctst)));
    bp(i)=mean(mean(imadjust(tf(:,:,3),ctst*2)));

%reference    
    a(i)=eval(['mean(mean(imadjust(rgb2gray(',eqb,':)),ctstr)));']);
    r(i)=eval(['mean(mean(imadjust((',eqb,'1)),ctstr+0.1)));']);
    g(i)=eval(['mean(mean(imadjust((',eqb,'2)),ctstr)));']);
    b(i)=eval(['mean(mean(imadjust((',eqb,'3)),ctstr+0.1)));']);

    i=i+1;
end
%%
figure;subplot(4,1,1);imshow(imadjust(uint8(frame(:,:,2,217*30)),ctst));subplot(4,1,2);image(uint8(frame(:,:,:,15*30)))
subplot(4,1,3);eval(['imshow(imadjust(',eqb,'2),ctstr));subplot(4,1,4);image(uint8(',eqb,':)))']);
% figure;imshow(uint8(-frame(:,:,3,15*30)+mean(frame(:,:,3,15*30),'all')));
%%
st=10*30;
x=1:length(r);
x=x/30;
figure;plot(x(st:end),mapminmax(ap(st:end)),x(st:end),mapminmax(a(st:end)));


% % %%
figure;
subplot(4,1,1);plot(x(st:end),zscore(r(st:end)),x(st:end),zscore(rp(st:end)));
subplot(4,1,2);plot(x(st:end),zscore(g(st:end)),x(st:end),zscore(gp(st:end)));
subplot(4,1,3);plot(x(st:end),zscore(b(st:end)),x(st:end),zscore(bp(st:end)));
subplot(4,1,4);plot(x(st:end),zscore(a(st:end)),x(st:end),zscore(ap(st:end)));

%%
starttime=find(zscore(g(st:end))<-1,1); %find the stimulation start time

span=0.5;hp=1;
[~,d]=lowpass(r,span,30,ImpulseResponse="iir",Steepness=0.95);
% first step normalization
[r1,rp1]=ampadj(r,rp,st,d);
[g1,gp1]=ampadj(g,gp,st,d);
[b1,bp1]=ampadj(b,bp,st,d);
[a1,ap1]=ampadj(a,ap,st,d);

% chose the second step normalization methods based on detection results
bt=g1(st:end)';
bpt=gp1(st:end)';



% bt=mapminmax(b1(st:end),0,1);
% bpt=mapminmax(bp1(st:end),0,1);
% bt=zscore(r1(st:end));
% bpt=zscore(rp1(st:end));
% bt=bt-mean(bt);
% bpt=bpt-mean(bpt);
t=1:length(bt);t=t/30;
figure;plot(t,bt,t,bpt);

% lowpass/ smooth for denoise
pb=bt-bpt;
% pb=g1(st:end)'-gp1(st:end)'+b1(st:end)'-bp1(st:end)'+r1(st:end)'-rp1(st:end)';
pb=smooth(pb,64,'sgolay');

[~,dpb]=lowpass(pb,0.33,30,ImpulseResponse="iir");
pb=filtfilt(dpb,pb);

figure;plot(t,smooth(bt,29),t,smooth(bpt,29),t,pb);

tP1=ones(starttime+411*30+st-1,1)*mean(pb);
tP1(st:st+length(pb)-1)=pb;
P1=tP1;
P3=P1(starttime+st-21*30-1-1:starttime+411*30+st-1);
Ps=P3;
t=1:length(Ps);t=t/30;
figure();plot(t,Ps);hold on; scatter(21,Ps(21));hold off;
save([date,enum,'_-21s.mat'],'Ps');
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

function [r1,rp1]=ampadj(r,rp,st,d)

r1 = filtfilt(d,r);rp1 = filtfilt(d,rp);

r =  filloutliers(r,"clip");
rp =  filloutliers(rp,"clip");
r1=r;rp1=rp;
% r1=smooth(r,10,'rlowess');rp1=smooth(rp,10,'rlowess');
end