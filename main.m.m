function main.m  % detection of copy-move forgery in images
clc
clear
clear all;
%% parameters
blocksize=8;
overlapp=1;
Nd=16;
Th=0.9999;    %from (0.9 - 1)
s_threshold=2;
%% input image
[name,path]=uigetfile('*.*');
img=imread([path,name]);
%img=imread('s.bmp');
[r c n] = size(img);
if n > 1
    im=rgb2gray(img);
else
    im=img; % It's already gray.
end
figure,subplot(1,3,1),imshow(im),title('Input Image')
%% divide into overlapping blocks
tic
a=1;
for j=1:overlapp:(c-blocksize)+1
    for i=1:overlapp:(r-blocksize)+1
        sondos(a).block=im(i:i+blocksize-1,j:j+blocksize-1);
        sondos(a).position=[i j];
        sondos(a).index=a;
        a=a+1;
    end
end
%% applay DCT for each block
sz=size(sondos,2);
DC=[];
QZ=4;
FDCT=[];
for a=1:sz
    [feature,vector]=featureExtraction(sondos(a).block);
    DC(a,1)=feature(1,1);
    FDCT(a,:)=feature;
end
%% divide into groups
G=[];
numclass=20;
[centers,mincenter,mindist,q2,quality] = FastKmean(FDCT,numclass,1);
for n=1:numclass
    ind= find(mincenter==n);
    for i=1:length(ind)
        G(n,i)= ind(i);
    end
end
%% draw segmentation
col=jet(size(G,1));
subplot(1,3,2),imshow(im),title('Clustered Image')
for e=1:size(G,1)
    color=col(e,:);
    for ee=1:size(G,2)
        if G(e,ee)~=0
            rectangle('position',[sondos(G(e,ee)).position(2),sondos(G(e,ee)).position(1) ,blocksize,blocksize],'Edgecolor',color)
        end
    end
end
%% detect copy move
subplot(1,3,3),imshow(im),title('Result Image')
for nG=1:size(G,1) % you can try parfor (parallel) instead of for
    emp=find( G(nG,:)==0);
    if isempty(emp)==0
        A=[];
        for a=1:emp(1)-1
            [f,vec]=featureExtraction(sondos(G(nG,a)).block);
            A(a,1:9)=f;
            A(a,10)=sondos(G(nG,a)).position(1);
            A(a,11)=sondos(G(nG,a)).position(2);
            
        end
    else
        A=[];
        for a=1:size(G,2)
            [f,vec]=featureExtraction(sondos(G(nG,a)).block);
            A(a,1:9)=f;
            A(a,10)=sondos(G(nG,a)).position(1);
            A(a,11)=sondos(G(nG,a)).position(2);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Asorted=RadixSort(A,9);
    for i=1:size(Asorted,1)-1
        similar=[];
        for l=1:9  %num of features
            s=abs(Asorted(i+1,l)-Asorted(i,l));
            if s<s_threshold
                %similar%
                similar(l)=1;
            else
                %not similar%
                similar(l)=0;
            end
        end
        if isempty(find(similar==0))    %two block is similar calculat distance
            x1=Asorted(i,10);
            x2=Asorted(i+1,10);
            y1=Asorted(i,11);
            y2=Asorted(i+1,11);
            D= sqrt((x1-x2)^2+(y1-y2)^2);
            if D>Nd  %calculate shift vector
                rectangle('Position',[y1,x1,blocksize,blocksize],'Edgecolor','r');
                rectangle('Position',[y2,x2,blocksize,blocksize],'Edgecolor','r');
                % line([y1,y2],[x1,x2],'Color','r','LineWidth',1)
            end
        end
    end  
end
time=toc;
end 