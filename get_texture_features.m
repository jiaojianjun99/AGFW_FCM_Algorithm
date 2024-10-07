function Texture = get_texture_features(I,label,N,G,Offset)
Texture=[];
for i=1:3
image=I(:,:,i);
% % I=rgb2gray(I);
% I0=I==0;
% I_coodinate=find(I0);%检测图像中是否有非零元素
 texture =[]; 
 for label_Num= 1:N  
     
    Region=label==label_Num;
    Region=uint8(Region);
      %转化成标签下的原始图像
     R = image.*Region;
     glcm=graycomatrixjiao(R, 'NumLevels', 8, 'G',[G(1),G(2)],'Offset', [Offset(1),Offset(2)]);%灰度共生矩阵  
%%    glcm =graycomatrixexplain( R, 'NumLevels', 8, 'G',[1,255],'Offset', [1  -1] );
     stats = graycoprops( glcm, 'Contrast Energy Homogeneity ' );   
     texture_v = [stats.Contrast,stats.Energy,stats.Homogeneity];
     texture = [ texture ; texture_v ];  
 end
 Texture=[Texture,texture];
 if i==2
     jiao=1;
 end
  
 
end
end


