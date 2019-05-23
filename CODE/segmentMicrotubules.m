function [segmentedData,microTubules_L,dataOut,dataOut2] = segmentMicrotubules(dataIn,distanceThreshold)
%function [segmentedData,microTubules_L,dataOut,dataOut2] = segmentMicrotubules(dataIn,distanceThreshold)
% This function is called from another function thus it will be assumed
% that the data is passed as a 3D matlab matrix. 

%% Parse input
if nargin<1
    help segmentMicrotubules;
    segmentedData=[];
    return
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
    help segmentMicrotubules;
    segmentedData=[];
    return
else
    % the data is not 4D, it has to be 3D as in [R,G,B].
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