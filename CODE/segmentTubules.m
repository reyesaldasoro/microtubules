function [cellTubules]=segmentTubules(dataIn,cellBody,cellProtrusions)
%%
[rows,cols,~,~]                        = size(dataIn);
% filter to obtain a slightly better segmentation
sizeFilter                              = 5;
filtG                                   = gaussF(sizeFilter,sizeFilter,1);

%% Find an initial segmentation of the cell bodies, they may be merged
currentData                             = dataIn(:,:,:,1);
channel_1                               = (currentData(:,:,1));          % select the first  level of the 3D matrix
channel_2                               = (currentData(:,:,2));          % select the second level of the 3D matrix
channel_2_noCells                       = channel_2.*imerode(uint8(1-cellBody),ones(3));
%channel_1F                              = imfilter((channel_1),filtG,'replicate');
%channel_2F                              = imfilter((channel_2),filtG,'replicate');
%%
% Otsu threshold
%level_2                                 = 255* graythresh(channel_2F);
%level_1                                 = 255* graythresh(channel_1F);


%% Only calculate the edges OUTSIDE the cells and close to them
%cellBody_3                              = imopen(cellBody_2,ones(3));
distFromCells                           = bwdist(cellBody+cellProtrusions);
regionCells                             = (distFromCells)<85;
regionEdges                             = regionCells-imdilate(cellBody+cellProtrusions,ones(3));

%imagesc(cellBody+2*cellNuclei)
%% Prepare for the tubules
% Detect the edges
BW                                      = edge(channel_2,'canny',[],0.75);
% remove all those that are far from the cells or on the cells themselves
BW2                                     = regionEdges.*(BW);
[BW2_L,numEdges]                        = bwlabel(BW2);

BW2_Props                               = regionprops(BW2_L,channel_2,'Area',...
    'MajoraxisLength','MinoraxisLength',...
    'Eccentricity','Euler','MaxIntensity','BoundingBox');


%% Iterate over all edges, find the pixels that surround them and take only that are brighter than the rest, 
% 50-50 for a line 33-66 for a horseshoe
BW3 = zeros(rows,cols);
for k=1:numEdges
    if (BW2_Props(k).Area>5)
        %disp(k)
        rr                  = max(1,-3+floor(BW2_Props(k).BoundingBox(2))):min(rows,3+ceil(BW2_Props(k).BoundingBox(2)+BW2_Props(k).BoundingBox(4)));
        cc                  = max(1,-3+floor(BW2_Props(k).BoundingBox(1))):min(cols,3+ceil(BW2_Props(k).BoundingBox(1)+BW2_Props(k).BoundingBox(3)));
        currentEdge         = bwmorph(BW2_L(rr,cc)==k,'thin');
        currentDist         = bwdist(currentEdge);
        currentBoundary     = (floor(currentDist)==1)+(currentDist==2);
        % important to remove the cell from the intensitites otherwise it can lead to
        % a very high threshold
        %currentIntens       = channel_2(rr,cc).*uint8(1-cellBody(rr,cc));
        currentIntens       = channel_2_noCells(rr,cc);
        currentHistog       = currentIntens(currentBoundary==1);
        currentThres        = 255*graythresh(currentHistog);
        currentRidge        = ((currentIntens>=(0.9*currentThres)).*(currentBoundary));
        currentRidge2       = (imdilate(currentRidge,ones(3)));
        currentRidge3       = bwmorph(currentRidge2,'thin','inf');

        %BW5(rr,cc) = BW5(rr,cc) + k*imfill(imclose(BW3(rr,cc)==k,ones(3)),'holes');
        BW3(rr,cc) = BW3(rr,cc) + k*(currentRidge3);
        %imagesc(BW5);drawnow;pause(0.05)
    end
end


% find endpoints of the edges to connect those that are separated as in 1 1
% 1 0 1 1 1
BW2_endp                                = (bwmorph(BW3>0,'endpoints'));
% Now bridge between the endpoints, this is better than closing or dilating
% as only connections between end points are formed.

%BW2_endp_b                              = (bwmorph(BW2_endp,'bridge'));
% Dilate with a cross and then apply a majority, single edges will stay the
% same, but those tha are close will become a H and will keep the bridge
BW2_endp_d                              = imdilate(BW2_endp,[0 1 0;1 1 1;0 1 0]);
BW2_endp_m                              = (bwmorph(BW2_endp_d,'majority'));

%% Thin the edges and remove spurious pixels
BW4 = bwmorph(((BW3>0)+BW2_endp_m)>0,'thin');

% Label to discard very small edges
[BW4_L,numEdges]                          = bwlabel(BW4);
BW4_P = regionprops(BW4_L,channel_2,'Area',...
    'MajoraxisLength','MinoraxisLength',...
    'Eccentricity','Euler','MaxIntensity','BoundingBox');
%imagesc(BW)
%%

cellTubules      = (ismember(BW4_L,find([BW4_P.Area]>5)));

% %%
% BW6 = (bwmorph(BW5,'thin','inf'));
% BW7 = bwlabel(bwmorph(BW6,'spur',0));
% 
% %BW6 = BW5 - BW2;
% %BW7 = bwlabel(BW6==1);
% 
% BW8 = regionprops(BW7,channel_2,'Area','MinIntensity','MeanIntensity','MaxIntensity',...
%     'MajoraxisLength','MinoraxisLength',...
%     'Eccentricity','Euler','BoundingBox');
% %% Final tubules
% % Brightest tubules, take the minimum levels and compare with the max of each tubule
% % and also discard anything less than 5 pixels
% brightTubLev     = ceil( mean(double([BW8.MinIntensity]))+0.5*std(double([BW8.MinIntensity])));
% cellTubules      = (ismember(BW7,find([BW8.MaxIntensity]>brightTubLev))).*(ismember(BW7,find([BW8.Area]>5)));

%
%brightTubLev    = 0.05*max([BW8.MaxIntensity]) +0.95*   mean([BW8.MaxIntensity]);
%BW9A             =(ismember(BW7,find([BW8.MaxIntensity]>brightTubLev))).*(ismember(BW7,find([BW8.Area]>15)));
% Long and straight
%BW9B             =(ismember(BW7,find([BW8.Eccentricity]>0.96))).*(ismember(BW7,find([BW8.Area]>6)));
%imagesc((BW7>0)+BW9B+2*BW9A)
%% Now include those that are 1-2 pixels away of the previously selected
%BW10             = bwdist(BW9A + BW9B);
%BW11            = unique(BW7.*(BW10<3));
%BW12            = ismember(BW7,BW11(2:end));
%BW13            = imclose(BW12,ones(1));
%imagesc((BW7>0)+(BW13))
%cellTubules         = BW13;
%cellTubules         = ismember(BW7,find([BW8.MajorAxisLength]>6));

