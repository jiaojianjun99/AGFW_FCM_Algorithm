%% This code was written by Jianjun Jiao in Lanzhou Jiaotong University.
% Email:jiaojianjun@lzjtu.edu.cn
% Last change: January, 2023
%%  
clear all;
clc;
close all;
%% ��������
cluster_n =3;
N=30;%
max_iter =100;      % Max. iteration 
G=[1,255];
Offset=[0 1;-1 1;-1 0;-1 -1];
%%
I =imread('image\wenli1.png');
%%  �����ػ������
[m,n,p]=size(I);
[label,N1] = superpixels(I,N);%N1Ϊʵ�ʵĳ����ظ��� 
BW = boundarymask(label);%������boundarymask�ú�����
% �����ؽ��
figure,imshow(imoverlay(I,BW,'cyan'),'InitialMagnification',100);  %cyan������ɫ��imshow(a,'InitialMagnification','fit')ͼƬ������
idx = label2idx(label);%$$$$�ҵ�����$$$$$$$$$$$��ȡ��ǩλ�ã�ԭͼ��������
%���㳬���ؾ�ֵ
outputImage = zeros(size(I),'like',I);%%% ���ɺ�Iһ����С��ά�ȵ������
output = zeros(N1,3);%%%%%����
sumo=zeros(N1,3);%NΪ����������sumoΪ�������ڲ�����ֵ�ĺ�1
for labelVal = 1:N1

    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+m*n;%��־�ڶ�ҳ���ݣ����۽磩
    blueIdx = idx{labelVal}+2*m*n;%��־����ҳ����
    sumo(labelVal,:)=[sum(I(redIdx)),sum(I(greenIdx)),sum(I(blueIdx))];%��������������������ֵ�ĺͣ����������������󣻷ŵ�һҳ
    output(labelVal,:) = [mean(I(redIdx)),mean(I(greenIdx)),mean(I(blueIdx))];%%%%%�����صľ�ֵ���ŵ�һҳ
    %%%����ֵ�������ͼ����ά��
    outputImage(redIdx) = mean(I(redIdx));
    outputImage(greenIdx) = mean(I(greenIdx));
    outputImage(blueIdx) = mean(I(blueIdx));
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,1), 'sheet1'); % ��dataд��test.xls�Ĺ�����sheet1��
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,2), 'sheet2');
% xlswrite('d:\matlab2.xlsx ', outputImage(:,:,3), 'sheet3');
end    
rgb(:,:,1)=output(:,1);
rgb(:,:,2)=output(:,2);
rgb(:,:,3)=output(:,3);
hsv = rgb2hsv(rgb);%%%rgb��hsv֮���ת��
hsv1=mat2gray(hsv(:,:,1));
hsv2=mat2gray(hsv(:,:,2));
hsv3=mat2gray(hsv(:,:,3));

% figure,imshow(outputImage,'InitialMagnification',80)%outputImage��ֵ��������صĽ��
%% ��ȡ������������
Texture = get_texture_features(I,label,N1,G,Offset);
Texture1 =[sum(Texture(:,[1,4,7]),2),sum(Texture(:,[2,5,8]),2),sum(Texture(:,[3,6,9]),2)]./3;
texture =mat2gray(Texture1);%���ݹ�һ��
%% ��ȡ����ͳ������
%%  FCM�ָ�
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
%%��ʾ�ָ���


%%�������صķָ���չ��ԭ����ͼ��
label3=[];
for labelVal=1:N1  %labelValΪ�����ر�ǩ1-N1
    id=idx{labelVal};%�õ������������ǩΪ1��N�������ض�Ӧԭʼͼ�����ص������
    label3(id)=label2(labelVal);%������㶼����Ϊһ�����࣬label2�ǳ����ط��࣬����������ԭʼͼ��������ȥ��label3Ϊ1*154401
end
r2 = reshape(label3, m,n);     % ����ת��ΪͼƬ��ʽ
figure,imshow(label2rgb(r2,'jet')) ;


% %%ͬһ�೬�������ֵ,�����ɫ��
% cluster=unique(label2);%ͳ�ƴ���
% number=histc(label2,cluster);%ͳ��ÿһ����ռ����
% S=[accumarray(label2,rgb(:,:,1)),accumarray(label2,rgb(:,:,2)),accumarray(label2,rgb(:,:,3))];
% %Z��ʾ���ģ���ƽ��ֵ
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
% SS=W(label3,:,:);%��ά����
% 
% smap= uint8((reshape(SS, m,n,p)));
% smap = imresize(smap, [320, 540]);
% set (gcf,'Position',[0,0,540,320]);  % ǰ�������崰������Ļ��λ�ã����������ڴ�С
% % ʹͼ������Ӧ��������
% figure,imshow(smap,'border','tight','initialmagnification','fit');
% % set (gcf,'Position',[0,0,640,320]);  % ǰ�������崰������Ļ��λ�ã����������ڴ�С
% 
% %%ʹͼ������Ӧ��������
% % figure,imshow(smap,'border','tight','initialmagnification','fit');
