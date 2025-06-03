%% This code was written by Jianjun Jiao in Lanzhou Jiaotong University.
% Email:jiaojianjun@lzjtu.edu.cn
% Last change: Oct, 2024
%This code is applicable to MATLAB 2018b
% Better segmentation results can be obtained by adjusting the superpixel parameter N and the cluster number cluster_n
%%  
clear all;
clc;
close all;
%% Parameter setting
cluster_n =3;
N=30;%
max_iter =100;      % Max. iteration 
G=[1,255];
Offset=[0 1;-1 1;-1 0;-1 -1];
%%
I =imread('wenli1.png');
%%  Superpixel region
[m,n,p]=size(I);
[label,N1] = superpixels(I,N);
BW = boundarymask(label);
% 
figure,imshow(imoverlay(I,BW,'cyan'),'InitialMagnification',100); 
idx = label2idx(label);
%Calculate the superpixel mean
outputImage = zeros(size(I),'like',I);
output = zeros(N1,3);
sumo=zeros(N1,3);
for labelVal = 1:N1

    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+m*n;
    blueIdx = idx{labelVal}+2*m*n;
    sumo(labelVal,:)=[sum(I(redIdx)),sum(I(greenIdx)),sum(I(blueIdx))];
    output(labelVal,:) = [mean(I(redIdx)),mean(I(greenIdx)),mean(I(blueIdx))];
    %%
    outputImage(redIdx) = mean(I(redIdx));
    outputImage(greenIdx) = mean(I(greenIdx));
    outputImage(blueIdx) = mean(I(blueIdx));
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,1), 'sheet1'); 
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,2), 'sheet2');
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,3), 'sheet3');
end    
rgb(:,:,1)=output(:,1);
rgb(:,:,2)=output(:,2);
rgb(:,:,3)=output(:,3);
hsv = rgb2hsv(rgb);
hsv1=mat2gray(hsv(:,:,1));
hsv2=mat2gray(hsv(:,:,2));
hsv3=mat2gray(hsv(:,:,3));


%% Get regional texture features
Texture = get_texture_features(I,label,N1,G,Offset);
Texture1 =[sum(Texture(:,[1,4,7]),2),sum(Texture(:,[2,5,8]),2),sum(Texture(:,[3,6,9]),2)]./3;
texture =mat2gray(Texture1);

%%  FCM
output=mat2gray(output);
data(:,:,1)= texture(:,1);
data(:,:,2)= texture(:,2);
data(:,:,3)= texture(:,3);
data(:,:,4)= output(:,1);
data(:,:,5)= output(:,2);
data(:,:,6)= output(:,3);
data(:,:,7)= hsv1;
data(:,:,8)= hsv2;
data(:,:,9)= hsv3;


data=[data(:,:,1),data(:,:,2),data(:,:,3),data(:,:,4),data(:,:,5),data(:,:,6),data(:,:,7),data(:,:,8),data(:,:,9)].*256;
% data=[data(:,:,1),data(:,:,2),data(:,:,3),data(:,:,4),data(:,:,5),data(:,:,6)].*256;
% data=[data(:,:,1),data(:,:,2),data(:,:,3)].*256;
% [~, label2]=STFRFCM(data,cluster_n,N1,3);
[~, label2,ss1]=AGFW_FCM(data,cluster_n,N1,max_iter);


%% The segmentation of superpixel is extended to the original image
label3=[];
for labelVal=1:N1  
    id=idx{labelVal};
    label3(id)=label2(labelVal);
end
r2 = reshape(label3, m,n);     
figure,imshow(label2rgb(r2,'jet')) ;


%%
cluster=unique(label2);
number=histc(label2,cluster);
S=[accumarray(label2,rgb(:,:,1)),accumarray(label2,rgb(:,:,2)),accumarray(label2,rgb(:,:,3))];
%
if cluster_n ==2;
Z(1,:)=S(1,:) /number(1);
Z(2,:)=S(2,:) /number(2);
end

if cluster_n ==3;
Z(1,:)=S(1,:) /number(1);
Z(2,:)=S(2,:) /number(2);
Z(3,:)=S(3,:) /number(3);
end

if cluster_n ==4;
Z(1,:)=S(1,:) /number(1);
Z(2,:)=S(2,:) /number(2);
Z(3,:)=S(3,:) /number(3);
Z(4,:)=S(4,:) /number(4);
end
if cluster_n ==5;
Z(1,:)=S(1,:) /number(1);
Z(2,:)=S(2,:) /number(2);
Z(3,:)=S(3,:) /number(3);
Z(4,:)=S(4,:) /number(4);
Z(5,:)=S(5,:) /number(5);
end
W(:,:,1)=Z(:,1);
W(:,:,2)=Z(:,2);
W(:,:,3)=Z(:,3);

SS=W(label3,:,:);%

smap= uint8((reshape(SS, m,n,p)));
smap = imresize(smap, [320, 540]);
set (gcf,'Position',[0,0,540,320]);  
% 
figure,imshow(smap,'border','tight','initialmagnification','fit');
% set (gcf,'Position',[0,0,640,320]); 


