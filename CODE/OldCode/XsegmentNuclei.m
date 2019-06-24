function nuclei = segmentNuclei(dataIn)
%%

channel_1                           = (dataIn(:,:,1));          % select the first  level of the 3D matrix
%%
level_1                             = 255* graythresh(channel_1);
% NUCLEI
% First do a thresholding, based on otsu and clean with a majority
nuclei_1                            = bwmorph(channel_1>(0.5*level_1),'majority');
nuclei_1P                           = regionprops(nuclei_1,'Area');

% nuclei can appear not completely solid, merge with a closing
nuclei_2                            = (imclose(nuclei_1,strel('disk',2)));
[nuclei_2L,numClasses]              = bwlabel(nuclei_2);

% discard very small dots (Area<10) and smaller (10<A<50) CLOSE to larger ones, but keep them if they are far away
nuclei_2P                           = regionprops(nuclei_2,'Area');
%
nuclei_3                            = ismember(nuclei_2L,find([nuclei_2P.Area]>10));
[nuclei_3L,numClasses]              = bwlabel(nuclei_3);

% discard very small dots (Area<10) and smaller (10<A<50) CLOSE to larger ones, but keep them if they are far away
nuclei_3P                           = regionprops(nuclei_3,'Area');


currentDistance(numClasses)         =0;

for k=1:numClasses
    %
    
    tempDistance                    = bwdist(nuclei_3L==k);
    distToNeighboursM               = (nuclei_3).*(tempDistance);
    distToNeighbours                = distToNeighboursM(nuclei_3>0);
    distToNeighbours(distToNeighbours==0)=[];
    distToNeighbours                = sort(distToNeighbours);
    currentDistance(1,k)              = distToNeighbours(1);
    
    %figure(2)
    %imagesc(tempDistance)
    
    %
end

smallNuclei                         = ismember(nuclei_3L,find([nuclei_3P.Area]<60));
closeNuclei                         = ismember(nuclei_3L,find(currentDistance(1,:)<50));
%currentDistance(2,:)=[nuclei_3P.Area];

nuclei                              = nuclei_3.*(1-smallNuclei.*closeNuclei);

% Finally do an erosion as low thresholds may prompt mergers

%imagesc(nuclei_2+nuclei)
%imagesc(nuclei)

%%

