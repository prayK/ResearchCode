clear
clc
row=1080;                                              %���ص����ó�1080pix
g1=zeros(row);                                         %����1080������󣬳�ʼ��g1
g2=zeros(row);                                         %����1080������󣬳�ʼ��g1
w0=100;                                                %�����뾶
l=3;                                                   %���˺�ֵ
for m1=1:row                                           %�����⺯��
for n1=1:row
if (m1-row/2).^2+(n1-row/2).^2<row/2*row/2
g1(n1,m1)=sqrt((m1-row/2)^2+(n1-row/2)^2)^abs(l)*exp(-((m1-row/2)^2+(n1-row/2)^2)/w0^2)*exp(1i*l*atan2((n1-row/2),(m1-row/2))); 
end
end
end
A0=max(max(abs(g1)));                                 %���沨����
for m1=1:row
for n1=1:row
if (m1-row/2).^2+(n1-row/2).^2<row/2*row/2
g2(n1,m1)=A0.*exp(-((m1-row/2)^2+(n1-row/2)^2)/(row*0.2)^2).*exp(1i*pi*-2*((m1-row/2)^2+(n1-row/2)^2)/(row/6)^2); 
end
end
end