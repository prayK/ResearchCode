clc
clear all
close all
%%
c=-916:916;
r=-916:916;
[x,y]=meshgrid(c,r);
[theta,r]=cart2pol(x,y);

figure;
for l=1:8
    subplot(2,4,l)
    g=mod(l*theta,2*pi);
    imshow(g,[])%自动调整数据的范围以便于显示
    l=l+1;
end

figure;
for l=1:8
    subplot(2,4,l)
    g=mod(-1*l*theta,2*pi);
    imshow(g,[])
    l=l+1;
end
%imshow(mod(-2*theta,2*pi),[])