clear all;close all;
Image_Input=imread('9.bmp');
dims=ndims(Image_Input);                     %��ȡ�����ά��
if(dims==3)
   Image_Gray=rgb2gray(Image_Input);           
elseif(dims~=2)
    return;
end
[h,w]=size(Image_Gray);                  
xr=8;
xl=8;
yt=8;
yb=8;
width=w-xr-xl;
height=h-yt-yb;
Image_Cut=imcrop(Image_Gray,[xr,yt,width-1,height-1]);   %  �г��ڱ�
figure(1);
imshow(Image_Cut,[]);
[Row,Column]=size(Image_Cut);
Image_FFT=fftshift(fft2(Image_Cut));
Image_Mod=abs(Image_FFT);
figure(1);
imshow(Image_Mod/Row/Column);
title('����ͼ��Ƶ��');
Max_Mod=max(max(Image_Mod));
Image_Ln=log10(Image_Mod/Max_Mod);  % ��һ����ȡ����
if rem(Row,2)~=0         %��ȡˮƽ���Ĵ�ƽƵ��
    Image_array=Image_Ln((Row+1)/2,1:end);      
else
    Image_array=(Image_Ln(Row/2,1:end)+Image_Ln(Row/2+1,1:end))/2;
end
if(rem(Column,2)==0)                            %ʹ���������Գ�������ʾ
    X=-Column/2:Column/2-1;
else
    X=-(Column-1)/2:(Column-1)/2;
end
figure(3);
plot(X,Image_array);
xlabel('��ɢƵ��')
ylabel('����ҶƵ�ʹ�һ����ȡ��Ȼ����');
if rem(Column,2)==0
    Data=Image_array(Column/2+1:end);
else
    Data=Image_array((Column+1)/2:end);
end
figure(4);
semilogx(Data,'b:*');
xlabel('��ɢƵ�ʵ���Ȼ����');
ylabel('����ҶƵ�ʹ�һ����ȡ��Ȼ����');
hold on;
if rem(Column,2)==0
    k=1:Column/2;
else
    k=1:Column/2+1;
end
p=polyfit(log10(k),Data,1);
line=polyval(p,log10(k));
semilogx(line,'r:*');
hold off;




 









        
