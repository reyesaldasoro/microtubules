%% Segmentation of Microtubules 

%% Clear all variables and close all figures
clear all
close all
clc 



%% Read the files that have been stored in the current folder
% These folders are specific for my computers (CCRA) so they need to be changed
% accordingly

if strcmp(filesep,'/')
    % Running in Mac    
    cd ('/Users/ccr22/OneDrive - City, University of London/Acad/Research/KCL/BrianStramer')
else
    % running in windows
    cd ('D:\OneDrive - City, University of London\Acad\Research\KCL\BrianStramer')
end

%% Read data and process with segmentMicrotubules
% If the data is stored as a multiple frame TIFF, they need to be converted to a 4-D
% matlab matrix, this can be done with the code below. 
%dataIn =readMTIFF('datasets/pair1.tif');
%dataIn =readMTIFF('datasets/2012_11_20_3.tif');

% If the data is stored with one tiff as red and one as green, then they are
% converted with this code
%[dataIn]=readMTIFF_G_R('/Users/ccr22/OneDrive - City, University of London/Acad/Research/KCL/BrianStramer/SingleCells/Control_ClipGFPRS_031018_3-1.tif');

% Otherwise, if they have already been saved as a 4-D matrix, just read
%load 2012_11_20_3
load('/Users/ccr22/OneDrive - City, University of London/Acad/Research/KCL/BrianStramer/SingleCells/Control_ClipGFPRS_031018_3-1.mat')
%load('/Users/ccr22/OneDrive - City, University of London/Acad/Research/KCL/BrianStramer/SingleCells/Control_ClipGFPRS_031018_1-2.mat')


    
%% Segment with the new process

% The number of timepoints is required later
numTimePoints       =size(dataIn,4);
% Select one time point
k=1;
% Each of these lines will perform a process and it is given by the name of the file.
% Times are calculated but these can be ignored or removed.

tic;[cellBody,cellNuclei,cellProtrusions,cellNoNuclei]  = segmentCellNuclei(dataIn(:,:,:,k));t1=toc;
tic;[clumps,notClumps,degreeClump,cellBody_L]           = analyseCellConditions(cellBody,cellNuclei);t2=toc;
tic;[cellTubules]                                       = segmentTubules(dataIn(:,:,:,k),cellBody,cellProtrusions);t3=toc;
tic;[cellTubules_L,cellBody_L_Complete]                 = allocateTubules(cellBody_L,cellProtrusions,cellTubules,cellNoNuclei);t4=toc;
tic;[dataOut_C, dataOut_CT,dataOut_CT2,dataOut_CT3]     = prepareDataOut(dataIn(:,:,:,k),cellBody_L_Complete,cellNuclei,cellTubules_L);t5=toc;
disp([t1 t2 t3 t4 t5 t1+t2+t3+t4+t5])
%% Display section 
% The results can be displayed in several ways: The original data with the cell,
% nuclei and tubules overlaid with colour code
imagesc(dataOut_CT3)

% Just the regions
% imagesc(cellBody+2*cellNuclei+ 4* cellProtrusions+ 5*cellTubules)

%% Increase contrast
% To appreciate the tubule segmentation, a higher contrast may be used
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
figure
% subplot(131)
% imagesc(dataIn(:,:,:,k))
% 
% subplot(132)
% imagesc(dataIn(:,:,2,k).*(1-uint8(imdilate( zerocross(cellNuclei-(cellBody+cellProtrusions)),ones(3)))))
% colormap(jet4)
% subplot(133)
%imagesc(cellNuclei+cellBody+0.3*cellProtrusions+0.5*cellTubules)
imagesc(double(dataIn(:,:,2,k)).*(1-(cellTubules>0)).*(1-(imdilate(zerocross(cellNuclei-(cellBody+cellProtrusions)),ones(3)))))
colormap(jet4)


%% Prepare to create a video
h0=gcf;
hDataOut = imagesc(dataOut_CT3);
hTitle   = title('');
clear F;
%% Iterate over the frames to create a video of the whole process
for k=1:numTimePoints
    
    [cellBody,cellNuclei,cellProtrusions,cellNoNuclei]  = segmentCellNuclei(dataIn(:,:,:,k));
    [clumps,notClumps,degreeClump,cellBody_L]           = analyseCellConditions(cellBody,cellNuclei);
    [cellTubules]                                       = segmentTubules(dataIn(:,:,:,k),cellBody,cellProtrusions);
    [cellTubules_L,cellBody_L_Complete]                 = allocateTubules(cellBody_L,cellProtrusions,cellTubules,cellNoNuclei);
    [dataOut_C, dataOut_CT,dataOut_CT2,dataOut_CT3]     = prepareDataOut(dataIn(:,:,:,k),cellBody_L_Complete,cellNuclei,cellTubules_L);
    hTitle.String                                       = strcat('Time =',32,num2str(k));
    hDataOut.CData                                      = dataOut_CT3;
    drawnow;
    pause(0.01)
    F(k) = getframe(h0);
end
%% Write the video
%v = VideoWriter('Tubules_2012_11_20_3.mp4', 'MPEG-4');
%v = VideoWriter('Control_ClipGFPRS_031018_1-2.mp4', 'MPEG-4');
v = VideoWriter('Control_ClipGFPRS_031018_3-1.mp4', 'MPEG-4');


open(v);
writeVideo(v,F);
close(v);
