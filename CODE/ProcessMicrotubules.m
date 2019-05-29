

%% Clear all variables and close all figures
clear all
close all
clc

%% Read the files that have been stored in the current folder
if strcmp(filesep,'/')
    % Running in Mac
    
    cd ('/Users/ccr22/OneDrive - City, University of London/Acad/Research/KCL/BrianStramer')
    %baseDir                             = 'Metrics_2019_04_25/metrics/';
else
    % running in windows
    %load('D:\OneDrive - City, University of London\Acad\ARC_Grant\Datasets\DataARC_Datasets_2019_05_03.mat')
    cd ('D:\OneDrive - City, University of London\Acad\Research\KCL\BrianStramer')
    %baseDir                             = 'Metrics_2019_04_25/metrics/';
end

%% Read data and process with segmentMicrotubules
%dataIn =readMTIFF('datasets/pair1.tif');
%dataIn =readMTIFF('datasets/2012_11_20_3.tif');
load 2012_11_20_3

%%
jet4 = [ 0         0    0
         0         0    0.6354
         0         0    0.7083
         0         0    0.7812
         0         0    0.8542
         0         0    0.9271
         0         0    1.0000
         0    0.2000    1.0000
         0    0.4000    1.0000
         0    0.6000    1.0000
         0    0.8000    1.0000
         0    1.0000    1.0000
    0.2500    1.0000    0.7500
    0.5000    1.0000    0.5000
    0.7500    1.0000    0.2500
    1.0000    1.0000         0
    1.0000    0.8333         0
    1.0000    0.6667         0
    1.0000    0.5000         0
    1.0000    0.3333         0
    1.0000    0.1667         0
    1.0000         0         0
    0.9881         0         0
    0.9762         0         0
    0.9643         0         0
    0.9524         0         0
    0.9405         0         0
    0.9286         0         0
    0.9167         0         0
    0.9048         0         0
    0.8929         0         0
    0.8810         0         0
    0.8690         0         0
    0.8571         0         0
    0.8452         0         0
    0.8333         0         0
    0.8214         0         0
    0.8095         0         0
    0.7976         0         0
    0.7857         0         0
    0.7738         0         0
    0.7619         0         0
    0.7500         0         0
    0.7381         0         0
    0.7262         0         0
    0.7143         0         0
    0.7024         0         0
    0.6905         0         0
    0.6786         0         0
    0.6667         0         0
    0.6548         0         0
    0.6429         0         0
    0.6310         0         0
    0.6190         0         0
    0.6071         0         0
    0.5952         0         0
    0.5833         0         0
    0.5714         0         0
    0.5595         0         0
    0.5476         0         0
    0.5357         0         0
    0.5238         0         0
    0.5119         0         0
    0.5000         0         0];


%%
tic;[cellBody,cellNuclei]               =segmentCellNuclei(dataIn);toc
%%
k=57;
[cellBody,cellNuclei,cellProtrusions]               =segmentCellNuclei(dataIn(:,:,:,k));
subplot(121)
imagesc(dataIn(:,:,:,k))

subplot(122)
imagesc(dataIn(:,:,2,k).*(1-uint8(imdilate( zerocross(cellNuclei-cellBody),ones(3)))))
colormap(jet4)



%%
tic;[cellBody,cellNuclei,cellTubules]               =segmentTubulesCellNuclei(dataIn);toc
tic;[segmentedData,microTubules,dataOut,dataOut2]   =segmentMicrotubules(dataIn);toc
%%

%figure
subplot(131)
imagesc(dataIn(:,:,2,1))
colormap(jet4)

subplot(132)
imagesc(dataIn(:,:,2,1).*(1-uint8(imdilate(cellTubules+ zerocross(cellBody),ones(2)))))
colormap(jet4)

subplot(133)
imagesc(dataIn(:,:,2,1).*(1-uint8(imdilate(microTubules,ones(2)))))
colormap(jet4)
%%
a=b;
[segmentedData,microTubules,dataOut,dataOut2]=segmentMicrotubules(a(:,:,:,68));
%%
kkk=38.416;
dataOut3 = kkk+(-kkk*repmat(nuclei,[1 1 3])+double(a(:,:,:,68)));

imagesc(dataOut3/max(dataOut3(:)))

%imagesc(a(:,:,:,68))
%%


[segmentedData,microTubules,dataOut,dataOut2]=segmentMicrotubules('datasets/pair1.tif');



%%


dir0 =('/Users/ccr22/Academic/work/microscopicCells/BrianStramer/datasets/pair1_mat_Tu/');
dir1 =dir(strcat(dir0,'*.mat'));

%%
load(strcat(dir0,dir1(1).name))
dataIn =(microTubules==0)+(microTubules==2);
imagesc(dataIn)
%%
 finalPointBits = calculatePointyEdges(dataIn);
 %imagesc(dataIn+finalPointBits)

 
 %%
 
load(strcat(dir0,dir1(50).name))
dataIn              =(microTubules==1)+(microTubules==2);
dataInSkel          = bwmorph(dataIn,'skel',inf);
dataInEroded        = imopen(imerode(dataIn,strel('disk',4)),ones(5));
pointyBits          = dataInSkel.*(dataIn-dataInEroded);

imagesc(dataIn+dataInEroded)

pointyBits_L        = bwlabel(pointyBits);
pointyBits_A        = regionprops(pointyBits_L,'Area','MajorAxisLength','MinorAxisLength');
largePointyB        = ismember(pointyBits_L,find([pointyBits_A.MajorAxisLength]>20));
largeThinPointyB    = ismember(pointyBits_L,find((([pointyBits_A.MajorAxisLength]<=20).*[pointyBits_A.MajorAxisLength]>12).*([pointyBits_A.MinorAxisLength]<5)));

imagesc(dataIn+dataInSkel+pointyBits+1*largePointyB+2*largeThinPointyB)
 
 

%% create the frames for a movie
for k=1:49
    
    %imagesc(dataOut2(:,:,:,k)/255)
    figure(1);imagesc(segmentedData(:,:,k))
    %figure(2);imagesc(microTubules(:,:,k))
    figure(3);imagesc(dataOut(:,:,:,k)/255)
    figure(4);imagesc(dataOut2(:,:,:,k)/255)
    drawnow;pause(0.05);
    %F(k) = getframe;
    
    
end


%% save the movie as avi
   movie2avi(F,'TubuleTracking3.avi','compression','none','fps',8)

    %% save the movie as a GIF
    [imGif,mapGif] = rgb2ind(F(1).cdata,256,'nodither');
    numFrames = size(F,2);

    imGif(1,1,1,numFrames) = 0;
    for k = 2:numFrames 
      imGif(:,:,1,k) = rgb2ind(F(k).cdata,mapGif,'nodither');
    end
   

    imwrite(imGif,mapGif,'TubuleTracking3.gif','DelayTime',0,'LoopCount',inf) %g443800












%%
tempa = a(:,:,:,4);


channel_1   = tempa(:,:,1);          % select the first  level of the 3D matrix
channel_2   = tempa(:,:,2);          % select the second level of the 3D matrix
channel_3   = tempa(:,:,3);          % select the third  level of the 3D matrix

imagesc(channel_2)
%%
level_1     =255* graythresh(channel_1);
level_2     =255* graythresh(channel_2);
level_3     =255* graythresh(channel_3);

nuclei      = channel_1>(1*level_1);
cellBody    = imerode(imclose(imopen((channel_2>(1*level_2)),strel('disk',1)),strel('disk',6)),strel('disk',0));
cellBody    = bwmorph(bwmorph(cellBody,'majority'),'majority');

imagesc(3*nuclei+2*cellBody.*(1-nuclei)+(channel_2>level_2))
%%
b(1,1,3) =0;



figure(1)
 subplot(211)
 imagesc(tempa)

c           = tempa;

c(:,:,1)    = c(:,:,1)+255*uint8(zerocross(-0.5+cellBody));
c(:,:,2)    = uint8(2*(double(c(:,:,2))))+255*uint8(zerocross(-0.5+nuclei));

c(c>255) = 255;
 subplot(212)
imagesc(c)
%%


level_2b     =255* graythresh(channel_2.*(uint8(cellBody)));


Tubules     = ((1-imdilate(cellBody,ones(5))).*double(channel_2)>0.78099*level_2b);
Tubules2     = bwmorph(Tubules,'clean');

Tubules3     = edge(channel_2,'canny',[],1.06623);
Tubules4     = (1-imdilate(cellBody,ones(1))).*Tubules3;

Tubules5    = bwmorph(Tubules4,'bridge');

%% Detect strong filaments with Canny, they have to be closed one 
% one side, both sides parallel, elongated and nothing else

%then detect very weak filaments, discard everything else and look for weak lines

%crossing filaments create a problem as they can generate squares of triangles


% figure(1)
% imagesc(Tubules4+Tubules5)

% Tubules5     = watershed(-imfilter(channel_2,gaussF(9,9,1),'replicate'));
% Tubules6     = (1-imdilate(cellBody,ones(5))).*Tubules5;

%imagesc(Tubules+Tubules2)
figure(2)

c(:,:,3)    = c(:,:,3)+255*uint8(zerocross(-0.5+Tubules5));

imagesc(double(channel_2)+25*(zerocross(-0.5+Tubules5)))
%%


cellBodyD           = imdilate(cellBody,strel('disk',5));
strongTubules       = ((res_strength>3.3).*(1-cellBodyD));
weakTubules         = imclose(bwlabel((res_strength>1.9965).*(1-cellBodyD)),ones(2));
weakNextStrong      = unique(strongTubules.*(weakTubules));

% keep weak tubules with a peak over the strong ones

finalTubules_1      = ismember(weakTubules,weakNextStrong(2:end));
finalTubules_2      = bwmorph(finalTubules_1,'thin','inf');
finalTubules_3      = bwmorph(finalTubules_2,'spur',1);
finalTubules_4      = bwmorph(finalTubules_3,'clean');

imagesc(double(channel_2)+50*finalTubules_4);
