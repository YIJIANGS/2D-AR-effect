clear;clc;close all
v = VideoReader('vid1.mp4');
[embedGIF,map]=imread('gif1.gif','frames','all'); 
% embed2=imresize(embed_im,[400 300]);
%%
% embed_im=ind2rgb(A(:,:,:,44),map);
% embed_im=embed_im(:,:,1);
% imshow(embed_im)
% for i=1:100
%     embed_im=ind2rgb(embedGIF(:,:,:,i),map);
% embed_im=embed_im(:,:,1);
% imshow(embed_im)
% pause(.1)
% end
%%
% v = VideoReader('vid1.mp4');
% while hasFrame(v)
%     im = readFrame(v);
%     imshow(im)
%     pause(.1)
% end
%%
% for i=1:2:800
%     embed_im=ind2rgb(embedGIF(:,:,:,i),map);
% embed_im=embed_im(:,:,1);
% imshow(embed_im)
% pause(.1)
% end
%%
embedFrame=50;
% while hasFrame(v)
for frame=1:50
    im = readFrame(v);
    embedFrame=embedFrame+4;
    embed_im=ind2rgb(embedGIF(:,:,:,embedFrame),map);
    embed_im=embed_im(:,:,1);
    embed2=imresize(embed_im',[400 300]);
%     g=rgb2gray(im);
bw=im(:,:,1)>135 & im(:,:,2)>110 & im(:,:,3)>100;
b2=bwareaopen(bw,13000);
b3=bwconvhull(b2);
[n,m]=find(b3);
if std(m)>std(n)
    [~,ul_idx]=min(n);
    [~,ur_idx]=max(m);
    [~,ll_idx]=min(m);
    [~,lr_idx]=max(n);
else
[~,ul_idx]=min(m+n);
[~,ur_idx]=max(m+size(b3,1)-n);
[~,ll_idx]=max(n+size(b3,2)-m);
[~,lr_idx]=max(n+m);
end
X=[m(ul_idx) m(ur_idx) m(lr_idx) m(ll_idx)];        %四个新顶点
Y=[n(ul_idx) n(ur_idx) n(lr_idx) n(ll_idx)];  

x=[100 400 400 100];        %四个原顶点（长方形）
y=[100 100 500 500];
% embed2=imresize(embed_im,[y(3)-y(1)+1,x(2)-x(1)+1]);
I2=zeros(size(b3));
I2(y(1):y(3)-1,x(1):x(2)-1)=embed2;
g2=I2<.5;
% I2=im2double(I2);

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';   %变换后的四个顶点，方程右边的值
%联立解方程组，方程的系数
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);             
  0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
  0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
  0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
  0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=inv(A)*B;        %用四点求得的方程的解，也是全局变换系数
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
m=fa(7);l=fa(8);

% [M, N] = size(bw);
newI=im2double(im);
k=0;
for i=x(1):x(3)
    for j=y(1):y(3)
        if g2(j,i)==1
            new_x=round((a*i+b*j+c)/(m*i+l*j+1));
            new_y=round((d*i+e*j+f)/(m*i+l*j+1));
            newI(new_y,new_x,:)=I2(j,i,:);
        end
    end
end

imshow(newI)
pause(.1)
[I,map2]=rgb2ind(newI,256);
    if frame==1
        imwrite(I,map2,'1.gif','gif', 'Loopcount',inf,'DelayTime',0.05);
    else
        imwrite(I,map2,'1.gif','gif','WriteMode','Append','DelayTime',0.05);
    end
    
end