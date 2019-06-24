function [segmentedData,microTubules_L,dataOut,dataOut2] = segmentMicrotubules(dataIn,dataInName,distanceThreshold)
%function [segmentedData,microTubules_L,dataOut,dataOut2] = segmentMicrotubules(dataIn,dataInName,distanceThreshold)

%% Parse input
if nargin<1
    help segmentMicrotubules;
    segmentedData=[];
    return
end

if isa(dataIn,'char')
    %data has been received as the name of the file to be processed,
    % Assume for now that it is an 4-dimensional TIFF File [rows,cols,levs,NUM_Timepoints]   
    dataInName                              = dataIn;
    dataIn                                  = readMTIFF(dataInName);
end

% Verify that there is a name
if ~exist('dataInName','var')
    dataInName                              = 'currentFileName.tif';
end

%verify that the distance has been set
if ~exist('distanceThreshold','var')
    distanceThreshold                       = 15;
end
%%
% find dimensions of the input data
[rows,cols,levs,timePoints]                 = size(dataIn);

%%
% if the data is 4D, then cycle over the timePoints otherwise process only
% the current time point
if timePoints>1
    dataIn4D                                = dataIn;
    segmentedData(rows,cols,timePoints)     = 0;
    microTubules_L(rows,cols,timePoints)    = 0;
    dataOut(rows,cols,levs,timePoints)      = 0;
    dataOut2(rows,cols,levs,timePoints)     = 0;
    
    % Name to be used to save the files in order with a number of zeros before the
    % number identifier
    dataOutName0                            = 'T0000';
    % Create Folders where the data will be saved as individualt time
    % points and in Matlab format
    dataFolder1                             = strcat(dataInName(1:end-4),'_mat_Or',filesep);
    dataFolder2                             = strcat(dataInName(1:end-4),'_mat_La',filesep);
    dataFolder3                             = strcat(dataInName(1:end-4),'_mat_Tu',filesep);
    [statF1,messF1,messIDF1]                = mkdir(dataFolder1);
    [statF2,messF2,messIDF2]                = mkdir(dataFolder2);
    [statF3,messF3,messIDF3]                = mkdir(dataFolder3);
    if ~isempty(messF1)
        %remove files from folder
        delete(strcat(dataFolder1,'*.*'));
    end
    if ~isempty(messF2)
        %remove files from folder
        delete(strcat(dataFolder2,'*.*'));
    end
    if ~isempty(messF3)
        %remove files from folder
        delete(strcat(dataFolder3,'*.*'));
    end
            
    for counterTime = 1:timePoints  
                        
        disp(strcat('Processing time point =',num2str(counterTime),'/',num2str(timePoints)));
        %[segmentedData(:,:,counterTime),microTubules_L(:,:,counterTime),dataOut(:,:,:,counterTime),dataOut2(:,:,:,counterTime)] = segmentMicrotubules(dataIn4D(:,:,:,counterTime),distanceThreshold);
        [microTubules,microTubules_L,dataOut,dataOut2] = segmentMicrotubules(dataIn4D(:,:,:,counterTime),dataInName,distanceThreshold);
        dataFileName                        = strcat(dataOutName0(1:end-floor(log10(counterTime))),num2str(counterTime));
        %dataOutName1            = strcat(dataOutFolder,dataOutName);
        %--------------------------------------------------------------------------
        %----- the images read are saved to a file HERE ---------------------------
        %----- the original data is saved as matlab
        dataIn                              = dataIn4D(:,:,:,counterTime);
        save(strcat(dataFolder1,dataFileName),'dataIn');
        %----- the Nuclei are saved as labelled data to be tracked separately
        dataL                               = bwlabel(microTubules==1);
        save(strcat(dataFolder2,dataFileName),'dataL');
        %----- all the segmentation: Cells, nuclei, tubules as a single image
        save(strcat(dataFolder3,dataFileName),'microTubules');
                                
    end
    
else
    % the data is not 4D, it can be 2D or 3D.
    if levs==3
        %% PERFORM SEGMENTATION OF THE NUCLEI AND THE CELLS 
        % WITH OTSU AND MORPHOLOGICAL OPERATORS
        channel_2                           = (dataIn(:,:,2));          % select the second level of the 3D matrix
     
        %cellNuclei                              = segmentNuclei(dataIn);
        %cellBody                            = segmentCellBody(dataIn,cellNuclei);
        [cellBody,cellNuclei]               = segment_Cell_Nuclei(dataIn);
        [cellNuclei_L,numNuclei]            = bwlabel(cellNuclei);
        [cellBody_L,numCells]               = bwlabel(cellBody);
%%        
        finalPointyBits                     = calculatePointyEdges(cellBody);
        microTubules_L                      = linearStructures(channel_2,cellBody);
        microTubules                        = (microTubules_L>0);
%%        
        [assignedTubules,finalCells,distFromCell] = assignTubulesToCells(microTubules_L,cellNuclei_L,cellBody_L);
        
        segmentedData                       = 1*cellNuclei+2*cellBody.*(1-cellNuclei)+3*assignedTubules.*(1-cellBody);
        % Test distance between tubules and OTHER cells
        % find the tubules that create the connection between cells and create output
        % D_CT          corresponds to the distance between cell and tubule
        % D_CC          corresponds to the distance between cells (that have been split)
        % D_Tubulues    corresponds to the tubules that have created the connection
        [D_CT,D_CC,D_TUBULES, dataOut,dataOut2]=calculateDistanceTubulesCell(numNuclei,assignedTubules,finalCells,distFromCell,microTubules_L,distanceThreshold,dataIn);
        
    else
        disp('Error: Data should be colour images with Nuclei=Red, Tubules=Green');
        return;
    end
end
%%
function [D_CT,D_CC,D_TUBULES,dataOut,dataOut2]=calculateDistanceTubulesCell(numNuclei,assignedNuclei,finalCells,distFromCell,microTubules_L,distanceThreshold,dataIn)

%%

D_CT                                        = zeros(numNuclei);
D_CC                                        = zeros(numNuclei);
D_TUBULES{numNuclei,numNuclei}              = 0;

for k1=1:numNuclei;
    %qq                                      = distFromCell(:,:,k1);
    currentDistances                        = distFromCell(:,:,k1);
    for k2=1:numNuclei;
        if k1~=k2
            
            
            distOverTubules                 = currentDistances.*(assignedNuclei==k2);
            connectingTubulesIm             = microTubules_L((distOverTubules<distanceThreshold)&(distOverTubules>0));
            connectingTubules               = unique(connectingTubulesIm);
            [yd1,xd1]                       = hist(currentDistances(assignedNuclei==k2),(1:200));
            [yd2,xd2]                       = hist(currentDistances(finalCells(:,:,k2)>0),(1:200));
            
            %plot(xd1,(yd1)/sum(yd1),'b-',xd2,(yd2)/sum(yd2),'r--')
            %grid on
            D_CT(k1,k2)                     = sum(yd1(1:distanceThreshold))/sum(yd1);
            D_CC(k1,k2)                     = sum(yd2(1:distanceThreshold))/sum(yd2);
            D_TUBULES{k1,k2}                = connectingTubules;
            
        end
    end
end
D_CT = D_CT .*(1-(D_CC>0));


    %% create output with tubules in colours according to nuclei
    
    dataOut                                     = dataIn;
    
        dataOut2                                    = double(dataIn);
    %% Detect direct contact between the cells
    [c1,c2]     = find(D_CC);
    if ~isempty(c1)
        connectingCells                         = sum(finalCells(:,:,c1),3);
        connectingCells                         = (imdilate(zerocross(connectingCells),ones(3)));
        
        dataOut2(:,:,3)                         = dataOut2(:,:,3)+(220*connectingCells);
        %imagesc(uint8(dataOut2))
    end
    %% Detect the tubule that creates the contact between the cells
    [c1,c2]     = find(D_CT);

    if ~isempty(c1)
        for counterContact = 1:numel(c1)
            %firstCell                                   = c1(counterContact);
            %secondCell                                  = c2(counterContact);
            contactTubules                              = D_TUBULES{c1(counterContact),c2(counterContact)}';
            currentTubules                              = ismember(microTubules_L,contactTubules);
            switch (rem(counterContact,7))
                case 0
                    tempdata                            = dataOut2(:,:,1);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,1)                     = tempdata;
                case 1
                    tempdata                            = dataOut2(:,:,2);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,2)                     = tempdata;
                case 2
                    tempdata                            = dataOut2(:,:,3);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,3)                     = tempdata;
                case 3
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                case 4
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,3);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,3)                     = tempdata2;
                    
                case 5
                    tempdata1                           = dataOut2(:,:,3);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,3)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                    
                case 6
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata3                           = dataOut2(:,:,3);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    tempdata3(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                    dataOut2(:,:,3)                     = tempdata3;
            end
            
            
            
            
        end
        
        
    end
    
    %imagesc(dataOut2/255)
    
    %% create output with tubules in colours according to nuclei
    
    
    for counterN=1:numNuclei
        switch (rem(counterN,7))
            case 0
                tempdata                            = dataOut(:,:,1);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata;
            case 1
                tempdata                            = dataOut(:,:,2);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,2)                      = tempdata;
            case 2
                tempdata                            = dataOut(:,:,3);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,3)                      = tempdata;
            case 3
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,2);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
            case 4
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,3);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,3)                      = tempdata2;
                
            case 5
                tempdata1                            = dataOut(:,:,3);
                tempdata2                            = dataOut(:,:,2);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,3)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
                
            case 6
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,2);
                tempdata3                            = dataOut(:,:,3);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
                dataOut(:,:,3)                      = tempdata3;
                
        end
        %     % assigns tubules divided in 3
        %     currentLevel                        = 1+rem(counterN,3);
        %     tempdata                            = dataOut(:,:,currentLevel);
        %     tempdata(assignedNuclei==counterN)  =255;
        %     dataOut(:,:,currentLevel)           = tempdata;
    end
    
    %imagesc(dataOut)
    dataOut2=uint8(dataOut2);

