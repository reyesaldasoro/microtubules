function [cellBody,cellNuclei,AllCells]= segment_Cell_Nuclei(dataIn)
% function [cellBody,cellNuclei]= segment_Cell_Nuclei(dataIn)
%
% segmentation process, first do an initial thresholding with otsu of the
% green channel, that will coarsely identify where the cells are, then use
% that as a mask to get regions where there can be nuclei, everything else,
% discard. 
% For each region identified, get the intensity of the nuclei, separately
% of others, that should be able to pick up very faint nuclei which are
% indeed located within a visible cell. Then do a segmentation of the
% nuclei. 
%% Find an initial segmentation of the cell bodies, they may be merged
channel_2                               = (dataIn(:,:,2));          % select the second level of the 3D matrix
[rows,cols]                             = size(channel_2);
% filter to obtain a slightly better segmentation
sizeFilter                              = 5;
channel_2F                              = imfilter((channel_2),gaussF(sizeFilter,sizeFilter,1),'replicate');
% Otsu threshold
level_2                                 = 255* graythresh(channel_2F);
% CELL BODY
cellBody_1F                             = channel_2F>(0.75*level_2);
cellBody_2                              = imfill(cellBody_1F,'holes');
cellBody_3                              = imerode(cellBody_2,strel('disk',3));
% imagesc(cellBody_3)
%% eliminate small sections that will not be cells, but rather tubules or noise
[cellBody_3L,numSegments]               = bwlabel(cellBody_3);
cellBody_3P                             = regionprops(cellBody_3L,'Area');

cellBody_sizes                          = sort([cellBody_3P.Area]);

% mean of the largest regions:
cellBody_means                          = mean(cellBody_sizes(ceil(0.75*numSegments):numSegments));

cellBody_4                              = ismember(cellBody_3L,find([cellBody_3P.Area]>(1+0.1*cellBody_means)));
[cellBody_4L,numLargeCells]             = bwlabel(cellBody_4);
cellBody_4P                             = regionprops(cellBody_4,'Area','BoundingBox');


% imagesc(cellBody_3+cellBody_4)
%% Search now for nuclei inside these regions
% This is good to detect some nuclei that can be dim and thus would be lost
% with a single threshold, but care is needed because the strong green
% fluorescence may leak into the red, so the level for the nuclei has to
% be considerably higher than surrounding areas to be considered as nuclei
channel_1                               = (dataIn(:,:,1));%.*uint8(cellBody_4));          % select the first  level of the 3D matrix

nuclei_1(rows,cols)                     = 0;
level_1(numLargeCells)                  = 0;
for counterSegments=1:numLargeCells
    rowsRegion                          = max(1,floor(cellBody_4P(counterSegments).BoundingBox(2))):min(rows,floor(cellBody_4P(counterSegments).BoundingBox(2))+cellBody_4P(counterSegments).BoundingBox(4));
    colsRegion                          = max(1,floor(cellBody_4P(counterSegments).BoundingBox(1))):min(cols,floor(cellBody_4P(counterSegments).BoundingBox(1))+cellBody_4P(counterSegments).BoundingBox(3));
    tempChannel                         = channel_1(rowsRegion,colsRegion);
    level_1(counterSegments)            = 255* graythresh(tempChannel);
    tempNuclei                          = (tempChannel>level_1(counterSegments));
    
    levNucl(counterSegments)            = mean(tempChannel(tempNuclei>0));
    levNoNuc(counterSegments)           = mean(tempChannel(tempNuclei==0));
    if (levNucl(counterSegments) >(10+levNoNuc(counterSegments) ))
    nuclei_1(rowsRegion,colsRegion)     = nuclei_1(rowsRegion,colsRegion) +tempNuclei;
    end
end


nuclei_1L                                = bwlabel(nuclei_1>0);
nuclei_1P                               = regionprops(nuclei_1L,'Area');

nuclei_2a                             = ismember(nuclei_1L,find([nuclei_1P.Area]>4));



% nuclei can appear not completely solid, merge with a closing

nuclei_2                                = (imclose(nuclei_2a,strel('disk',1)));
%%nuclei_2                                = bwmorph(nuclei_2a,'majority');


% discard very small dots (Area<m-3std) of the largest nuclei regions  IFF they are close to other nuclei

[nuclei_2L,numNuclei]                   = bwlabel(nuclei_2);
nuclei_2P                               = regionprops(nuclei_2L,'Area','BoundingBox');

nuclei_sizesNS                          = ([nuclei_2P.Area]);


%%
sizeDisk                                = 17;
%dilatingElement                         = strel('disk',sizeDisk);
dilatingElement                         = ones(sizeDisk);
vectorNuclei                            = 1:numNuclei;
discardNuclei                           = zeros(1,numNuclei);
for counterNuc=1:numNuclei
%counterNuc=9;
    rowsRegion                          = max(1,-sizeDisk+floor(nuclei_2P(counterNuc).BoundingBox(2))):min(rows,sizeDisk+floor(nuclei_2P(counterNuc).BoundingBox(2))+nuclei_2P(counterNuc).BoundingBox(4));
    colsRegion                          = max(1,-sizeDisk+floor(nuclei_2P(counterNuc).BoundingBox(1))):min(cols,sizeDisk+floor(nuclei_2P(counterNuc).BoundingBox(1))+nuclei_2P(counterNuc).BoundingBox(3));

    reducedNuclei                       = nuclei_2L(rowsRegion,colsRegion);
%    currentNuc                          = imdilate(nuclei_2L==counterNuc,dilatingElement);
%    otherNuclei                         = ismember(nuclei_2L,setdiff(vectorNuclei,counterNuc));
    currentNuc                          = imdilate(reducedNuclei==counterNuc,dilatingElement);
    otherNuclei                         = ismember(reducedNuclei,setdiff(vectorNuclei,counterNuc));
    combinationNuc                      = (otherNuclei+2*currentNuc);
    %imagesc(combinationNuc)
%    presentNuclei                       = unique(nuclei_2L(currentNuc>0));
    presentNuclei                       = unique(reducedNuclei(currentNuc>0));
    presentNuclei                       = presentNuclei(2:end);
    if numel(presentNuclei>1)
        % There are two nuclei close to each other, discard one IFF it is
        % much smaller than the other
        presentNuclei_Sizes             =  nuclei_sizesNS(presentNuclei);
        largestNuclei                   = 0.2*max(presentNuclei_Sizes);
        smallestNuclei                  = min(presentNuclei_Sizes);
        if (largestNuclei>smallestNuclei)
            discardNuclei(presentNuclei(presentNuclei_Sizes<largestNuclei))=1;
            
        end
    end 
end
nuclei_3                                = ismember(nuclei_2L,find(1-discardNuclei));

%% Keep the option to remove things near the edge (where the edge may be inside the image)

nuclei_4                                = (cellBody_4L>0).*(imclose(nuclei_3,strel('disk',2)));


% imagesc(nuclei_3+cellBody_4)



%imagesc(nuclei_3+nuclei_4)

%%

%% Find a segmentation of the cell bodies that may be merged

%simple and crude would be to use the watershed of the nuclei, but this may
%not follow the true boundaries of the cells.  Besides, close cells that
%may be open in a "C" shape.

[cellNuclei_L,numNuclei]                    = bwlabel(nuclei_4);
%[cellBody_L,numCells]                   = bwlabel(cellBody);

%detect which cells have more than one nuclei
%%
sizeDisk                                = 9;
dilatingElement                         = strel('disk',sizeDisk-2);
%%
cellBody_5(rows,cols)                   = 0;                      
for counterCell=1:numLargeCells
%%    %traverse each cell and decide over it
%counterCell=15;
    % Act only on the actual cell and nuclei of the cell
    rowsRegion                          = max(1,-sizeDisk+floor(cellBody_4P(counterCell).BoundingBox(2))):min(rows,sizeDisk+floor(cellBody_4P(counterCell).BoundingBox(2))+cellBody_4P(counterCell).BoundingBox(4));
    colsRegion                          = max(1,-sizeDisk+floor(cellBody_4P(counterCell).BoundingBox(1))):min(cols,sizeDisk+floor(cellBody_4P(counterCell).BoundingBox(1))+cellBody_4P(counterCell).BoundingBox(3));
     
    reducedCell                         = cellBody_4L(rowsRegion,colsRegion)==counterCell;
    reducedNuclei                       = cellNuclei_L(rowsRegion,colsRegion).*reducedCell; %remove other nuclei
    %figure(1);     imagesc(reducedCell+reducedNuclei)
    
    % Detect if there are more than one nuclei in the cell
    presentNuclei                       = unique(reducedNuclei(reducedCell));
    presentNuclei                       = presentNuclei(2:end);
    if numel(presentNuclei)==1
        % just one cell, for cells that do not overlap, do a closing of the cell
        cellBody_5(rowsRegion,colsRegion) = cellBody_5(rowsRegion,colsRegion) +imclose(reducedCell,dilatingElement);
        %figure(2);    imagesc(reducedCell+reducedNuclei+imclose(reducedCell,dilatingElement))       
    elseif numel(presentNuclei)==0
        % Cell has no nuclei, DISCARD?
    else
       % Split cells that overlap
       % Simplest option, split using only the nuclei, a Voronoi separation
       % or a watershed, But probably a better option is to fit 
       % separate Self-Organising Maps to each nuclei and then let to self
       % adapt
       
       reducedIntensity                     = imdilate(reducedCell,ones(1)).*double(channel_2(rowsRegion,colsRegion,:));
       % Use this to segment the cell with S.O.M.
       %[reducedCellSplit]                  = splitCells_SOM(reducedCell,reducedNuclei,reducedIntensity,level_2);
       
       %Use this to segment with a watershed based on the nuclei use a
       %watershed, then make sure that only 2 objects are returned and not
       %small parts
       splitRegions                             = (watershed(-reducedNuclei));
       boundaryCells                            = imdilate((splitRegions==0),[0 1 0;1 1 1;0 1 0]);
       finalSegmentation                        = (1-boundaryCells).*reducedCell;
       numRedNuclei                             = numel(unique(reducedNuclei(:)))-1;
       [finalSegmentation_L,numRedSegments]     = bwlabel(finalSegmentation>0);
        if numRedSegments>(numRedNuclei)
            finalSegmentation_P                 = regionprops(finalSegmentation_L,'Area');
            [largestSegments,indexSegments]     = sort([finalSegmentation_P.Area]);
            finalSegmentation                   = ismember(finalSegmentation_L,indexSegments(end-numRedNuclei+1:end));
        % beep;
        end

       
       
       cellBody_5(rowsRegion,colsRegion)    = cellBody_5(rowsRegion,colsRegion) + finalSegmentation;
       
       
    end
%%    
    
    
end
% imagesc(cellBody_3)
%%
cellBody                                = cellBody_5;
cellNuclei                              = nuclei_4;
AllCells                                = cellBody_2;
