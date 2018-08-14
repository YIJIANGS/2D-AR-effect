clear;clc;close all
im=imread('3.jpg');
embed_im=imread('i2.jpg');
imshow(im)
%%
g=rgb2gray(im);
bw=g>150;
imshow(bw)
%%
% b2=imdilate(bw,ones(5));
b2=bwareaopen(bw,10000);
b3=bwconvhull(b2);
imshow(b3)
%%
[n,m]=find(b3);
[~,ul_idx]=min(m+n);
[~,ur_idx]=max(m+size(b3,1)-n);
[~,ll_idx]=max(n+size(b3,2)-m);
[~,lr_idx]=max(n+m);
imshow(b3)
hold on
scatter(m(ul_idx),n(ul_idx))
scatter(m(ur_idx),n(ur_idx))
scatter(m(ll_idx),n(ll_idx))
scatter(m(lr_idx),n(lr_idx))
%%
X=[m(ul_idx) m(ur_idx) m(lr_idx) m(ll_idx)];        %四个新顶点
Y=[n(ul_idx) n(ur_idx) n(lr_idx) n(ll_idx)];  

x=[m(ul_idx) m(ur_idx) m(ur_idx) m(ul_idx)];        %四个原顶点（长方形）
y=[n(ul_idx) n(ul_idx) n(lr_idx) n(lr_idx)];
%%
embed2=imresize(embed_im,[y(3)-y(1)+1,x(2)-x(1)+1]);
imshow(embed2)
%%
I2=im;
I2(y(1):y(3),x(1):x(2),:)=embed2;
I2=im2double(I2);
%%

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
            new_x=round((a*i+b*j+c)/(m*i+l*j+1));
            new_y=round((d*i+e*j+f)/(m*i+l*j+1));
            newI(new_y,new_x,:)=I2(j,i,:);
    end
end
imshow(newI)
%%
% i=897;j=692;
% new_x=round((a*i+b*j+c)/(m*i+l*j+1))
% new_y=round((d*i+e*j+f)/(m*i+l*j+1))
%%
b3=zeros(size(bw));
b3(y(1):y(3),x(1):x(3))=1;
[m1,n1]=find(b3);
%%
new_x=round((a*m1+b*n1+c)./(m*m1+l*n1+1));
new_y=round((d*m1+e*n1+f)./(m*m1+l*n1+1));
new_ind=sub2ind(size(bw),new_y,new_x);
for i=1:3
    tempI=newI(:,:,i);
    tempI2=I2(y(1):y(3),x(1):x(3),i);
    tempI(new_ind)=0;
    newI(:,:,i)=tempI;
end
imshow(newI)