function cellBody =segmentCellBody(dataIn,nuclei)

%%

% use the nuclei to find boundaries
boundaryNuclei_1                        = watershed(-nuclei);

channel_2                               = (dataIn(:,:,2));          % select the second level of the 3D matrix
% filter to obtain a slightly better segmentation
sizeFilter                              = 5;
channel_2F                              = imfilter((channel_2),gaussF(sizeFilter,sizeFilter,1),'replicate');
% Otsu threshold
level_2                                 = 255* graythresh(channel_2F);
% CELL BODY
cellBody_1F                             = channel_2F>(0.75*level_2);
cellBody_2                              = imerode(cellBody_1F,strel('disk',1));

cellBody_2                              = imclose(cellBody_2,strel('disk',4));

cellBody_3                              = imfill(cellBody_2,'holes');
cellBody_3                              = imerode(cellBody_3,strel('disk',1));

% remove cells without nuclei
cellBody_3L                             = bwlabel(cellBody_3);
cellBody_3L_wNuclei                     = unique(cellBody_3L(nuclei>0));
cellBody_3L_wNuclei(cellBody_3L_wNuclei==0)=[];
cellBody_4                              = ismember(cellBody_3L,cellBody_3L_wNuclei).*(boundaryNuclei_1>0);

%%
cellBody_5                              = imerode(cellBody_4,strel('disk',1));
%cellBody                            = imfill(cellBody_4,'holes');
%imagesc(2*cellBody_2+cellBody_1)
%cellBody_3F                          = bwmorph(cellBody_2F,'majority');


%imagesc(cellBody_1F*2+cellBody_5)
%



%
cellBody_L                          = bwlabel(bwmorph(bwmorph(cellBody_5,'majority'),'majority'));
cellBody_size                       = regionprops(cellBody_L,'Area');
cellBody                            = ismember(cellBody_L,find([cellBody_size.Area]>90));
% VERIFY THAT THERE ARE NO CELLS THAT MERGE, IF SO SPLIT

%imagesc(cellBody+(boundaryNuclei_1==0))

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % enlargedNuclei                      = imdilate(nuclei,strel('disk',14));
% % sizeFilter=13;
% % CellIntensityFiltered                = imfilter(double(channel_2),gaussF(sizeFilter,sizeFilter,1));
% % boundaryNuclei_2                        = watershed(-CellIntensityFiltered.*(enlargedNuclei));
% % imagesc(2*(boundaryNuclei_1==0)+(enlargedNuclei)+(boundaryNuclei_2==0).*(enlargedNuclei))
% 
% %%
% % boundaryIntensity_1                     = watershed((-CellIntensityFiltered).*(1-nuclei)-100*nuclei);
% % imagesc(cellBody_3F+2*nuclei+2*(boundaryNuclei_1==0)+4*(boundaryIntensity_1==0))
% 
% 
% 
% %imagesc((boundariesIntensity==0)+2*(boundaryNuclei_1==0).*cellBody_3F)
% %%
%     %currentCellIntensity2           = (imfilter((40+double(channel_2)).*(1-currentCell),gaussF(sizeFilter,sizeFilter,1)));
%     %currentCellIntensity3           = currentCellIntensity2+currentCellIntensity;
% boundariesIntensity                 = watershed(-CellIntensityFiltered);
% 
% imagesc((boundariesIntensity==0)+2*cellBody_3F)
% 
% 
% %%
% 
% 
% boundaryNuclei_low                      = double(channel_2).*(watershedNuclei==0);
% boundaryNuclei_low(boundaryNuclei_low>(0.66*level_2)) =0;
% boundaryNuclei_high                     = double(channel_2).*(watershedNuclei==0);
% boundaryNuclei_high(boundaryNuclei_high<(0.66*level_2)) =0;
% boundaryNuclei_low                      =boundaryNuclei_low>0;
% boundaryNuclei_high                     =boundaryNuclei_high>0;
% 
% 
% imagesc((cellBody_3.*(boundaryNuclei_1>0))+2*(boundaryNuclei_high)+4*boundaryNuclei_low)
% %%
% %imagesc(cellBody_3+3*nuclei+(boundariesIntensity==0)+2*(boundaryNuclei_low))
% %
% newBoundaries     =  (boundariesIntensity==0).*bwdist(boundaryNuclei_low);
% %imagesc((newBoundaries<22).*(newBoundaries>0)+2*(boundaryNuclei_high))
% imagesc(nuclei+2*(boundaryNuclei_1==0)+(boundariesIntensity==0).*(imdilate(nuclei,strel('disk',21))))
% 
% %%
% cellBody_25                         = imfill(cellBody_2,'holes');
% cellBody_3                          = imclose(cellBody_25,strel('disk',5));
% cellBody_4                          = imerode(cellBody_3,strel('disk',1));
% cellBody                            = imfill(cellBody_4,'holes');
% imagesc(2*cellBody_2+cellBody_1)
% %%
% cellBody_L                          = bwlabel(bwmorph(bwmorph(cellBody,'majority'),'majority'));
% cellBody_size                       = regionprops(cellBody_L,'Area');
% cellBody                            = ismember(cellBody_L,find([cellBody_size.Area]>50));
% % VERIFY THAT THERE ARE NO CELLS THAT MERGE, IF SO SPLIT
% % the number is not enough as there can be cells in the edge that have no nuclei
% 
% 
% 
% % if unique((nuclei_L>0).*cellBody_L) is not sequential, there is acell with no
% % nuclei
% % how to detect a cell with 2 nuclei immediately????
% % duplicates in ((nuclei_L>0).*cellBody_L)
% 
% 
% 
% imagesc(  double(channel_2).*(cellBody))

