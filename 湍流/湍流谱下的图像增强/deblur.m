clear all;close all;
Image_Input=imread('19.bmp');
dims=ndims(Image_Input);                     %获取矩阵的维数
if dims==3
   Image_Gray=rgb2gray(Image_Input);           
elseif dims~=2
    return;
end
[h,w]=size(Image_Gray);                  
xr=8;
xl=8;
yt=8;
yb=8;
width=w-xr-xl;
height=h-yt-yb;
Image_Cut=imcrop(Image_Gray,[xr,yt,width-1,height-1]);   %  切除黑边
figure(1);
imshow(Image_Cut,[]);
[Row,Column]=size(Image_Cut);
Image_fft=fftshift(fft2(Image_Cut));    %傅里叶变换,平移
Image_abs=abs(Image_fft);    %求频谱
Image_norm=Image_abs/max(max(Image_abs));     %频谱归一化
Image_ln=log10(Image_norm);     %对数缩放处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image_array=Image_ln(fix(Row/2),1:end);      %获取水平中心处平频谱
% if rem(Column,2)==0                            %使用中心零点对称坐标显示
%     X=-fix(Column/2)+1:1:fix(Column/2);
% else
%     X=-fix(Column/2):1:fix(Column/2);
% end
% figure(3);
% plot(X,Image_array);
% xlabel('离散频率')
% ylabel('傅里叶频率归一化后取自然对数');
% Datablur=Image_array(fix(Column/2):end);
% figure(4);
% semilogx(Datablur,'b:*');
% xlabel('离散频率的自然对数');
% ylabel('傅里叶频率归一化后取自然对数');
% hold on;
% X_half=fix(Column/2):Column;
% count=size(X_half);
% k=1:count(2);
% p=polyfit(log10(k),Datablur,1);
% line=polyval(p,log10(k));
% semilogx(line,'r:*');
% hold off;
if rem(Row,2)~=0         %获取水平中心处平频谱
    Image_array=Image_ln((Row+1)/2,1:end);      
else
    Image_array=(Image_ln(Row/2,1:end)+Image_ln(Row/2+1,1:end))/2;
end
if rem(Column,2)==0                            %使用中心零点对称坐标显示
    Data_X=-Column/2:Column/2-1;
else
    Data_X=-(Column-1)/2:(Column-1)/2;
end
figure(2);
plot(Data_X,Image_array);
xlabel('离散频率')
ylabel('傅里叶频率归一化后取自然对数');
if rem(Column,2)==0
    Data_Y=Image_array(Column/2+1:end);
else
    Data_Y=Image_array((Column+1)/2:end);
end
figure(3);
semilogx(Data_Y,'b:*');
xlabel('离散频率的自然对数');
ylabel('傅里叶频率归一化后取自然对数');
hold on;
if rem(Column,2)==0
    k=1:Column/2;
else
    k=1:(Column+1)/2;
end
p=polyfit(log10(k),Data_Y,1);
line=polyval(p,log10(k));
semilogx(line,'r:*');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 直接将双对数坐标下的拟合直线斜率设为1  重建清晰图像的频谱  并计算参数
slope=-1;      % 给定双对数最佳拟合曲线的斜率，为大量统计得到的经验值
focal_point=7;  % 表示清晰图像频谱与退化图像频谱在双对数坐标系的下焦点力原点的距离
% 假设清晰图像频谱在双对数坐标系下的直线方程为 Y=aX+b,a是斜率，b是截距（待求解）
% Y=ln(G(0,focal_point)),X=ln(foacl_point),b=Y-aX
a=slope;
if rem(Column,2)==0
    x_length=Column/2;
    Y=Image_array(Column/2+focal_point);
    NewSpectrumHalf=zeros(1,Column/2);
else x_length=(Column-1)/2;
     Y=Image_array((Column+1)/2+focal_point);
     NewSpectrumHalf=zeros(1,(Column-1)/2);
end
X=log10(focal_point);
b=Y-slope*X;
for i=1:x_length
    x_temp=log10(i);
    if i<focal_point
        if rem(Column,2)==0
            NewSpectrumHalf(i)=Image_array(Column/2+i); 
        else NewSpectrumHalf(i)=Image_array((Column+1)/2+i);
        end
    else
        NewSpectrumHalf(i)=a*x_temp+b;
    end
end
New_spectrum_1=zeros(1,Column);
for i=1:Column                           
    if rem(Column,2)==0
        if i<=Column/2
            New_spectrum_1(i)=NewSpectrumHalf(1+Column/2-i);
        else
            New_spectrum_1(i)=NewSpectrumHalf(i-Column/2);
        end
    else
        if i<(Column+1)/2
            New_spectrum_1(i)=NewSpectrumHalf((1+Column)/2-i);
        elseif i==(Column+1)/2
            New_spectrum_1(i)=Image_array(i);
        else
            New_spectrum_1(i)=NewSpectrumHalf(i-(Column+1)/2);
        end
    end
end
figure(4);   
plot(Data_X,Image_array);
hold on;
plot(Data_X,New_spectrum_1,'r:*');
hold off;
Dvalue_1=Image_array-New_spectrum_1;
figure(5);
plot(Data_X,Dvalue_1);
xlabel('离散频率'); 
ylabel('含参数的代数式');
scale=100;
Data_X_stract=[fix(Column/2)-scale:fix(Column/2)+scale]-fix(Column/2);
Dvalue_stract_1=Dvalue_1(Data_X_stract+fix(Column/2));
figure(6);
plot(Data_X_stract,Dvalue_stract_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 用直线重建清晰图像的频谱y=ax+b  并计算参数
Number1=10;                                                  
Number2=20;                                               
Number3=5;                                                   % Number3用于修正直线的斜率
y1temp=Image_array(fix(Column/2)-Number1:fix(Column/2));     
y_1=sum(y1temp)/(Number1+1);                                 % 注意这里的数据个数不要搞错了，是Number1+1个
x1=fix(Column/2)-Number3;
x11=fix(Column/2)+1+Number3;
y2temp=Image_array(1:Number2);
y_2=sum(y2temp)/Number2+0.6;                                 % 若采用单张图像进行估计需要将末端整体上台一个1   因为退化前后末端的频谱不可能还是一样。                            
x2=1;                                                        
x22=Column;
coeff_a=(y_1-y_2)/(x1-x2);                     % 对称轴左边
coeff_b=((y_1+y_2)-coeff_a*(x1+x2))/2;
coeff_c=(y_1-y_2)/(x11-x22);                   % 对称轴右边
coeff_d=((y_1+y_2)-coeff_a*(x11+x22))/2;
for i=1:Column
    if i<=x1
         New_spectrum_2(i)=coeff_a*i+coeff_b;
    elseif i>x1 && i<x11
         New_spectrum_2(i)=Image_array(i);
    elseif i>=x11
         New_spectrum_2(i)=coeff_a*(Column-i)+coeff_b;
    end
end
figure(7);
plot(Data_X,Image_array);
hold on;
plot(Data_X,New_spectrum_2,'r*');
xlabel('离散频率');                             
ylabel('图像频谱归一化的自然对数');
hold off
Dvalue_2=Image_array-New_spectrum_2;
figure(8);
plot(Data_X,Dvalue_2);
xlabel('离散频率'); 
ylabel('含参数的代数式');
Dvalue_stract_2=Dvalue_2(Data_X_stract+fix(Column/2));
figure(9);
plot(Data_X_stract,Dvalue_stract_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha=0.001;
alpha=0.0008195 ;
belta=5/6;
Image_deblur=ones(Row,Column);
if rem(Row,2)==0
    for u=-Row/2:Row/2-1
        if rem(Column,2)==0
            for v=-Column/2:Column/2-1
                H_Estimate(u+Row/2+1,v+Column/2+1)=exp(-alpha*(u^2+v^2)^belta);
            end
        else 
            for v=-(Column-1)/2:(Column-1)/2
                H_Estimate(u+Row/2+1,v+(Column-1)/2+1)=exp(-alpha*(u^2+v^2)^belta);
            end
        end
    end
else
    for u=-(Row-1)/2:(Row-1)/2
        if rem(Column,2)==0
            for v=-Column/2:Column/2-1
                H_Estimate(u+(Row-1)/2+1,v+Column/2+1)=exp(-alpha*(u^2+v^2)^belta);
            end
        else
            for v=-(Column-1)/2:(Column-1)/2
                H_Estimate(u+(Row-1)/2+1,v+(Column-1)/2+1)=exp(-alpha*(u^2+v^2)^belta);
            end
        end
    end
end
% Image_deblur=Degeneration./H_Estimate;
% 
% Image_deblur=ifft2(ifftshift(Image_deblur));
% Image_deblur=real(Image_deblur);
% Image_deblur=0.5*(abs(Image_deblur)+Image_deblur);
% figure(11);
% imshow(Image_deblur,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_estimate=real(ifft2(ifftshift(H_Estimate))); % 转化到空域上来
result=deconvwnr(Image_Cut,h_estimate,0.001);
result=ifftshift(result); % 对图像进行1、3象限对调，2、4象限对调
figure(10);
imshow(result,[]);
title('维纳滤波复原图像');




