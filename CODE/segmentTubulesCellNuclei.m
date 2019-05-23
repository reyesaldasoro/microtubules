function [cellBody,cellNuclei,cellTubules]=segmentTubulesCellNuclei(dataIn)

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


% CELL BODY
cellBody_1F                             = channel_2F>(0.75*level_2);
cellBody_2                              = imfill(cellBody_1F,'holes');
cellBody_3                              = imopen(cellBody_2,ones(6));
cellBody_4                              = bwdist(cellBody_3);
cellBody                              = (cellBody_4)<60;

% Nuclei
cellNuclei                              = (channel_1F>(0.75*level_1));

%%

BW                                      = edge((uint8(1-cellBody_3).*channel_2),'canny',[],1);
BW2                                     = cellBody.*(BW.*imerode((1-cellBody_3),ones(3)));
[BW3,numEdges]                          = bwlabel(BW2);
BW4 = regionprops(BW3,channel_2,'Area',...
    'MajoraxisLength','MinoraxisLength',...
    'Eccentricity','Euler','MaxIntensity','BoundingBox');

%%
BW5 = zeros(rows,cols);
%%
for k=1:numEdges
%k=17;
    if (BW4(k).Area<180)&&(BW4(k).Area>10)
        %disp(k)
        rr = max(1,floor(BW4(k).BoundingBox(2))):min(rows,ceil(BW4(k).BoundingBox(2)+BW4(k).BoundingBox(4)));
        cc = max(1,floor(BW4(k).BoundingBox(1))):min(cols,ceil(BW4(k).BoundingBox(1)+BW4(k).BoundingBox(3)));
        
        BW5(rr,cc) = BW5(rr,cc) + imclose(BW3(rr,cc)==k,ones(3));
           
    end
 imagesc(BW5)
end
%%
BW6 = BW5 - BW2;
BW7 = bwlabel(BW6==1);

BW8 = regionprops(BW7,channel_2,'Area',...
    'MajoraxisLength','MinoraxisLength',...
    'Eccentricity','Euler','MaxIntensity','BoundingBox');
%% Final tubules
% Brightest tubules
brightTubLev    = 0.05*max([BW8.MaxIntensity]) +0.95*   mean([BW8.MaxIntensity]);
BW9A             =(ismember(BW7,find([BW8.MaxIntensity]>brightTubLev))).*(ismember(BW7,find([BW8.Area]>15)));
% Long and straight
BW9B             =(ismember(BW7,find([BW8.Eccentricity]>0.96))).*(ismember(BW7,find([BW8.Area]>6)));
%imagesc((BW7>0)+BW9B+2*BW9A)
%% Now include those that are 1-2 pixels away of the previously selected
BW10             = bwdist(BW9A + BW9B);
BW11            = unique(BW7.*(BW10<3));
BW12            = ismember(BW7,BW11(2:end));
BW13            = imclose(BW12,ones(3));

%imagesc((BW7>0)+(BW13))
cellTubules         = BW13;

