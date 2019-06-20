

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

%% Segment with the new process

k=3;
%tic;%[cellBody,cellNuclei]               =segmentCellNuclei(dataIn);toc
tic;[cellBody,cellNuclei,cellProtrusions]       = segmentCellNuclei(dataIn(:,:,:,k));t1=toc;
tic;[clumps,notClumps,degreeClump,cellBody_L]   = analyseCellConditions(cellBody,cellNuclei);t2=toc;
tic;[cellTubules]                               = segmentTubules(dataIn(:,:,:,k),cellBody,cellProtrusions);t3=toc;
disp([t1 t2 t3])
%toc
 imagesc(cellBody+2*cellNuclei+ 4* cellProtrusions+ 5*cellTubules)

%% display new process
figure
subplot(131)
imagesc(dataIn(:,:,:,k))

subplot(132)
imagesc(dataIn(:,:,2,k).*(1-uint8(imdilate( zerocross(cellNuclei-(cellBody+cellProtrusions)),ones(3)))))
colormap(jet4)
subplot(133)
%imagesc(cellNuclei+cellBody+0.3*cellProtrusions+0.5*cellTubules)
imagesc(double(dataIn(:,:,2,k)).*(1-(cellTubules>0)).*(1-(imdilate(zerocross(cellNuclei-(cellBody+cellProtrusions)),ones(3)))))
colormap(jet4)

