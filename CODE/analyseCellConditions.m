function [clumps,notClumps,degreeClump,cellBody_L]=analyseCellConditions(cellBody,cellNuclei)
%%

% Determine cells with 2 or more nuclei (clumps)
% First find where the nuclei are located
cellBody_L                              = bwlabel(cellBody);
nucleiPresent_centroids                 = (bwmorph(cellNuclei,'shrink','inf'));
nucleiPresent                           = cellBody_L(nucleiPresent_centroids);

numNuclei                               = numel(nucleiPresent);
indexNuclei                             = 1:numNuclei;
% Now determine where are more than 1 nuclei
% Find Unique ones
[C,ia,~]                               = unique(nucleiPresent);
% remove the positions of each of the uniques from the index
indexNuclei(ia)                         = [];
clumps                                  = unique(nucleiPresent(indexNuclei));
notClumps                               = setdiff(C,clumps);
%% Solid Clumps and just touching
% Further, detect those that are very close and only in contact by a very
% thin connection, find the degree of the connection
% The erosion may break a section of the cell but stay as a clump
if ~isempty(clumps)
numClumps = numel(clumps);
 
degreeClump(numClumps) = 0;
for counterClump = 1:numClumps
    numCellClump                    = 1;
    maxLabelNuclei                  = 1;
    sizeStrel                       = 1;
    %    while numCellClump==1
    while (maxLabelNuclei<2)&&(sizeStrel<35)
        sizeStrel                   = sizeStrel+1;
        [res,numCellClump]          = bwlabel(imerode(cellBody_L==clumps(counterClump),strel('disk',sizeStrel)));
        %labelNuclei                 = res(nucleiPresent_centroids);
        labelNuclei                 = unique(res.*cellNuclei);
        %maxLabelNuclei              = max(labelNuclei);
        maxLabelNuclei              = -1+numel(labelNuclei);
        degreeClump(counterClump)   = sizeStrel;
        %imagesc(res)
        %drawnow
    end
    %    imagesc(res)
end
else
    degreeClump =[];
    
end



