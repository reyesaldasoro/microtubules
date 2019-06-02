function [cellBody,cellNuclei,cellProtrusions]=segmentCellNuclei(dataIn)
%%
[rows,cols,levs,timeFrames]                        = size(dataIn);
% filter to obtain a slightly better segmentation
sizeFilter                              = 5;
filtG                                   = gaussF(sizeFilter,sizeFilter,1);
%% Find an initial segmentation of the cell bodies, they may be merged
currentData                             = dataIn(:,:,:,1);
channel_1                               = (currentData(:,:,1));          % select the first  level of the 3D matrix
channel_2                               = (currentData(:,:,2));          % select the second level of the 3D matrix
channel_1F                              = imfilter((channel_1),filtG,'replicate');
channel_2F                              = imfilter((channel_2),filtG,'replicate');
%%
% Otsu threshold
level_2                                 = 255* graythresh(channel_2F);
level_1                                 = 255* graythresh(channel_1F);


% Nuclei
cellNuclei_1                            = (channel_1F>(0.75*level_1));
cellNuclei                              = imclose(cellNuclei_1,strel('disk',3));

% CELL BODY
% Segment by intensity
cellBody_1                             = channel_2F>(0.95*level_2);
% Fill in the holes of the cells
cellBody_2                              = imfill(cellBody_1| imdilate(cellNuclei,strel('disk',3)),'holes');

% Some cells can be rather open, it is important to close the gaps of those
% cells, but only if: 1) they have only one nucleus (avoid clumps and cells
% without a nucleus) 2) they are within a range of sizes >400 pixels to
% avoid small specks 

cellBody_2L                             = bwlabel(cellBody_2);
% Remove small specks
cellBody_2P                             = regionprops(cellBody_2L,'Area');
cellBody_3                              = ismember(cellBody_2L,find([cellBody_2P.Area]>400));
cellBody_3L                             = bwlabel(cellBody_3);

% Remove cells without nuclei
cellsBody_3L_Nuclei                     = cellBody_3L.*cellNuclei;
cellsWithNuclei                         = unique(cellsBody_3L_Nuclei);
cellBody_4                              = ismember(cellBody_3L,cellsWithNuclei(2:end));
cellBody_4L                             = bwlabel(cellBody_4);

% Determine cells with 2 nuclei (clumps)
nucleiPresent_centroids                 =(bwmorph(cellNuclei,'shrink','inf'));
nucleiPresent_location                  = find(nucleiPresent_centroids);
nucleiPresent                           = cellBody_4L(nucleiPresent_location);



%%
cellBody_4E                             = imopen(cellBody_4L,strel('disk',3));
cellBody_4EP                            = regionprops(cellBody_4E,'Area','ConvexHull','ConvexImage','BoundingBox','Solidity','Eccentricity');
cellBody_4P                             = regionprops(cellBody_4L,'Area','ConvexHull','ConvexImage','BoundingBox','Solidity','Eccentricity');
% The combination for a hollow cell is that 1) it has a solidity below 0.8
% 2) it is not too large nor too small ( 700-2000) is not too elongated
% eccentricity < 0.9 and only one nucleus
%imagesc(cellBody_4E)

%%
cellBody_4R = cellBody_4L;%zeros(rows,cols);

for k=1:numel(cellBody_4EP)
    %k=17;
    if (cellBody_4EP(k).Solidity<0.8)
        if (cellBody_4EP(k).Area<2500)
            if (cellBody_4EP(k).Area>500)
                if (cellBody_4EP(k).Eccentricity<0.95)
                    if (sum(nucleiPresent==k)==1)
                        %disp(k)
                        rr = max(1,ceil(cellBody_4EP(k).BoundingBox(2))):min(rows,floor(cellBody_4EP(k).BoundingBox(2)+cellBody_4EP(k).BoundingBox(4)));
                        cc = max(1,ceil(cellBody_4EP(k).BoundingBox(1))):min(cols,floor(cellBody_4EP(k).BoundingBox(1)+cellBody_4EP(k).BoundingBox(3)));
                        
                        cellBody_4R(rr,cc) = (cellBody_4L(rr,cc)) + cellBody_4EP(k).ConvexImage;
                    end
                end
            end
        else
                                    rr = max(1,ceil(cellBody_4EP(k).BoundingBox(2))):min(rows,floor(cellBody_4EP(k).BoundingBox(2)+cellBody_4EP(k).BoundingBox(4)));
                        cc = max(1,ceil(cellBody_4EP(k).BoundingBox(1))):min(cols,floor(cellBody_4EP(k).BoundingBox(1)+cellBody_4EP(k).BoundingBox(3)));


            cellBody_4R(rr,cc) = (cellBody_4L(rr,cc)) + imfill(imclose(cellBody_4L(rr,cc),strel('disk',(5))),'holes');
        end
    end
end
cellBody_5                              = cellBody_4R>0;
%% Remove stems that are themselves tubules
cellBody_6                             = bwlabel(imopen(cellBody_5,strel('disk',5)));
cellBody_6P                             = regionprops(cellBody_6,'Area','ConvexHull','ConvexImage','BoundingBox','Solidity','Eccentricity');

cellBody_7                              = ismember(cellBody_6,find([cellBody_6P.Area]>400));

%%


cellBody                                = cellBody_7;
cellProtrusions_1                         = bwlabel(cellBody_5-cellBody_7);
cellProtrusions_P                            = regionprops(cellProtrusions_1,'Area','ConvexHull','ConvexImage','BoundingBox','Solidity','Eccentricity','majoraxislength');
cellProtrusions                              = ismember(cellProtrusions_1,find([cellProtrusions_P.MajorAxisLength]>10));
%imagesc(cellBody+2*cellNuclei)
