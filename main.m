%% This code was written by Jianjun Jiao in Lanzhou Jiaotong University.
% Email:jiaojianjun@lzjtu.edu.cn
% Last change: January, 2023
%%  
clear all;
clc;
close all;
%% 参数设置
cluster_n =3;
N=30;%
max_iter =100;      % Max. iteration 
G=[1,255];
Offset=[0 1;-1 1;-1 0;-1 -1];
%%
I =imread('image\wenli1.png');
%%  超像素获得区域
[m,n,p]=size(I);
[label,N1] = superpixels(I,N);%N1为实际的超像素个数 
BW = boundarymask(label);%￥￥￥boundarymask好函数！
% 超像素结果
figure,imshow(imoverlay(I,BW,'cyan'),'InitialMagnification',100);  %cyan，蓝绿色；imshow(a,'InitialMagnification','fit')图片的缩放
idx = label2idx(label);%$$$$我的手在$$$$$$$$$$$获取标签位置（原图以列排序）
%计算超像素均值
outputImage = zeros(size(I),'like',I);%%% 生成和I一样大小，维度的零矩阵
output = zeros(N1,3);%%%%%作用
sumo=zeros(N1,3);%N为超像素数，sumo为求超像素内部像素值的和1
for labelVal = 1:N1

    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+m*n;%标志第二页数据（开眼界）
    blueIdx = idx{labelVal}+2*m*n;%标志第三页数据
    sumo(labelVal,:)=[sum(I(redIdx)),sum(I(greenIdx)),sum(I(blueIdx))];%超像素区域内所有像素值的和；以列向量索引矩阵；放到一页
    output(labelVal,:) = [mean(I(redIdx)),mean(I(greenIdx)),mean(I(blueIdx))];%%%%%超像素的均值；放到一页
    %%%将均值赋给输出图像（三维）
    outputImage(redIdx) = mean(I(redIdx));
    outputImage(greenIdx) = mean(I(greenIdx));
    outputImage(blueIdx) = mean(I(blueIdx));
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,1), 'sheet1'); % 将data写入test.xls的工作表sheet1中
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,2), 'sheet2');
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,3), 'sheet3');
end    
rgb(:,:,1)=output(:,1);
rgb(:,:,2)=output(:,2);
rgb(:,:,3)=output(:,3);
hsv = rgb2hsv(rgb);%%%rgb到hsv之间的转换
hsv1=mat2gray(hsv(:,:,1));
hsv2=mat2gray(hsv(:,:,2));
hsv3=mat2gray(hsv(:,:,3));

% figure,imshow(outputImage,'InitialMagnification',80)%outputImage均值替代超像素的结果
%% 获取区域纹理特征
Texture = get_texture_features(I,label,N1,G,Offset);
Texture1 =[sum(Texture(:,[1,4,7]),2),sum(Texture(:,[2,5,8]),2),sum(Texture(:,[3,6,9]),2)]./3;
texture =mat2gray(Texture1);%数据归一化
%% 获取区域统计特征
%%  FCM分割
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
%%显示分割结果


%%将超像素的分割拓展到原来的图像
label3=[];
for labelVal=1:N1  %labelVal为超像素标签1-N1
    id=idx{labelVal};%得到的是这个（标签为1或N）超像素对应原始图像像素的坐标点
    label3(id)=label2(labelVal);%该坐标点都设置为一个分类，label2是超像素分类，将分类贴到原始图像像素上去，label3为1*154401
end
r2 = reshape(label3, m,n);     % 反向转化为图片形式
figure,imshow(label2rgb(r2,'jet')) ;


% %%同一类超像素求均值,最后颜色用
% cluster=unique(label2);%统计簇数
% number=histc(label2,cluster);%统计每一簇所占数量
% S=[accumarray(label2,rgb(:,:,1)),accumarray(label2,rgb(:,:,2)),accumarray(label2,rgb(:,:,3))];
% %Z表示中心，即平均值
% if cluster_n ==2;
% Z(1,:)=S(1,:) /number(1);
% Z(2,:)=S(2,:) /number(2);
% end
% 
% if cluster_n ==3;
% Z(1,:)=S(1,:) /number(1);
% Z(2,:)=S(2,:) /number(2);
% Z(3,:)=S(3,:) /number(3);
% end
% 
% if cluster_n ==4;
% Z(1,:)=S(1,:) /number(1);
% Z(2,:)=S(2,:) /number(2);
% Z(3,:)=S(3,:) /number(3);
% Z(4,:)=S(4,:) /number(4);
% end
% if cluster_n ==5;
% Z(1,:)=S(1,:) /number(1);
% Z(2,:)=S(2,:) /number(2);
% Z(3,:)=S(3,:) /number(3);
% Z(4,:)=S(4,:) /number(4);
% Z(5,:)=S(5,:) /number(5);
% end
% W(:,:,1)=Z(:,1);
% W(:,:,2)=Z(:,2);
% W(:,:,3)=Z(:,3);
% 
% SS=W(label3,:,:);%三维索引
% 
% smap= uint8((reshape(SS, m,n,p)));
% smap = imresize(smap, [320, 540]);
% set (gcf,'Position',[0,0,540,320]);  % 前两个定义窗口在屏幕的位置，后两个窗口大小
% % 使图像自适应填满窗口
% figure,imshow(smap,'border','tight','initialmagnification','fit');
% % set (gcf,'Position',[0,0,640,320]);  % 前两个定义窗口在屏幕的位置，后两个窗口大小
% 
% %%使图像自适应填满窗口
% % figure,imshow(smap,'border','tight','initialmagnification','fit');

