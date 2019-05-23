function finalPointyBits = calculatePointyEdges(cellBody)

%%
% Pad with zeros to avoid problems of cells in the boundary of the image
dataIn                                  = padData(cellBody,2,[],0);
%imagesc(paddedCells)
%%
% Label to process individually later
[dataIn_L, numCells]                    = bwlabel((dataIn));

%%
% Detect very large and pointy edges
sizeToErode                             = 7;
dataInThin                              = bwmorph(dataIn,'thin',sizeToErode);
dataInEroded                            = (imerode(dataIn,strel('disk',sizeToErode-2)));
dataInBoundary                          = dataIn-imerode(dataIn,[0 1 0;1 1 1;0 1 0]);
%
VlargeSize                              = sizeToErode*2;
minimumSize                             = 5;
allPointyBits_L                         = bwlabel((1-dataInEroded).*dataInThin);
allPointyBits_P                         = regionprops(allPointyBits_L,'Area');       
VlongPointyB                            = ismember(allPointyBits_L,find([allPointyBits_P.Area]>VlargeSize));
smallPointyB                            = ismember(allPointyBits_L,find(([allPointyBits_P.Area]<=VlargeSize).*([allPointyBits_P.Area]>minimumSize)));


imagesc(dataIn+smallPointyB+2*VlongPointyB)
%%
% Detect the pointy bits that lie within the boundary, then shrink them to
% a point, then these will be used to detect the angles at those points of
% the boundary, those with angles small will still classify as pointy bits
[pointyBitsToAnalyse,numPoints]         = bwlabel(bwmorph(dataInBoundary.*smallPointyB,'shrink','inf'));
%imagesc(dataIn+(pointyBitsToAnalyse))

%%
% loop over all pointy bits to Analyse

[rPoint,cPoint]                             = find(pointyBitsToAnalyse);
contourCells                                = bwboundaries(dataIn_L);
anglesPointyBits(numPoints,1)               = 0;
deltaV                                      = 1e-10;
for counterPoints = 1:numPoints
    %
    for distanceFromPoint                   = 7:20
        % Detect in which cell the point is located
        currentPoint                        = [rPoint(counterPoints),cPoint(counterPoints)];
        currentCell                         = dataIn_L(currentPoint(1),currentPoint(2));
        % extract contour for the current cell
        currentContour                      = contourCells{currentCell};
        numPointsBoundary                   = size(currentContour,1);
        % find the position of the pointy bit within the contour
        [q1,currPosition]                   = min(sum(abs(currentContour-repmat(currentPoint,[numPointsBoundary 1])),2));
        % shift the contour so that pointy bit is in between the other
        % points [1 pointybitPosition  2*pointybitPosition] to use to
        % calculate the slopes
        currentContour2                      =circshift(currentContour,-currPosition+distanceFromPoint+1);
        % calculate slopes of the two lines
        slope1                              = (currentContour2(1,2)-cPoint(counterPoints))/(deltaV+currentContour2(1,1)-rPoint(counterPoints));
        slope2                              = (currentContour2(2*distanceFromPoint+1,2)-cPoint(counterPoints))/(deltaV+currentContour2(2*distanceFromPoint+1,1)-rPoint(counterPoints));
        % tan(angle) = (m1 - m2)/(1+m1m2)
        tanAngle                            = -(slope1-slope2)/(1+slope1*slope2);
        currentAngle                        = 180*atan(tanAngle)/pi;
        % compensate for angles>90
        if (currentAngle<0)
            currentAngle                    = 180+currentAngle;
        end
        %disp([counterPoints currentAngle])
        %set(hpoint,'xdata',currentPoint(2),'ydata',currentPoint(1));
        % Store for each distance from the point
        anglesPointyBits(counterPoints,distanceFromPoint)     = round(currentAngle);
    end
end
%

%%
% The best way to estimate the angle is to get the MAXIMUM angle, as it may
% be small for short/long distances but we need ones that are consistent
% along the variation of the distance.
avAnglesPointyBits          = mean(anglesPointyBits(:,7:end),2);
anglesPointyBits(:,1)       = mean(anglesPointyBits(:,7:end),2);
anglesPointyBits(:,2)       = min(anglesPointyBits(:,7:end),[],2);
anglesPointyBits(:,3)       = max(anglesPointyBits(:,7:end),[],2);




% Only those angles <89 will be considered as pointy.
pointyBitsToKeep            = ismember(pointyBitsToAnalyse,find((anglesPointyBits(:,3))<89)); 

pointyBitsToKeep2           =unique(allPointyBits_L.*pointyBitsToKeep);

smallBitsToKeep             = ismember(allPointyBits_L,pointyBitsToKeep2(2:end));
imagesc(smallBitsToKeep + dataIn+2*VlongPointyB) 
%%
clf
imagesc(dataIn+2*imdilate(pointyBitsToAnalyse>0,ones(5)))
hold on
for counterPoints = 1:numPoints
    if anglesPointyBits(counterPoints,3)<60
        
        
        plot(cPoint(counterPoints),rPoint(counterPoints),'or','markersize',10)
        
    elseif anglesPointyBits(counterPoints,3)<80
        plot(cPoint(counterPoints),rPoint(counterPoints),'wd','markersize',10)
        
    elseif anglesPointyBits(counterPoints,3)<89
        plot(cPoint(counterPoints),rPoint(counterPoints),'gs','markersize',10)
        
    end
end


%%



% % Remove small regions 
% areaCells                             = regionprops(dataIn_L,'Area');
% dataIn                                = (ismember(dataIn_L,find([areaCells.Area]>250)));

%[dataIn_L, numCells]    = bwlabel(dataIn(:));
 


%paddedCells         = padData((microTubules==1)+(microTubules==2),3,[],0);
% Pad with zeros to avoid problems of cells in the boundary of the image
%paddedCells             = padData(dataIn,0,[],0);


% Obtain Skeletons, when pointy, a branch of the skeleton should traverse the point
dataInSkel              = bwmorph(dataIn,'skel',inf);

% erode the data so that only the edges of the skeleton remain
dataInEroded            = imopen(imerode(dataIn,strel('disk',4)),ones(5));
pointyBits              = dataInSkel.*(dataIn-dataInEroded);

% of the points that remain, label to find the axes lengths, 
% look a) for very long and b) long and thin
pointyBits_L            = bwlabel(pointyBits);
pointyBits_A            = regionprops(pointyBits_L,'MajorAxisLength','MinorAxisLength');
largePointyB            = ismember(pointyBits_L,find([pointyBits_A.MajorAxisLength]>20));
largeThinPointyB        = ismember(pointyBits_L,find((([pointyBits_A.MajorAxisLength]<=20).*[pointyBits_A.MajorAxisLength]>12).*([pointyBits_A.MinorAxisLength]<5)));

finalPointyBits         = largePointyB+largeThinPointyB;
%imagesc(dataIn+dataInSkel+1*largePointyB+1*largeThinPointyB)


% PREVIOUS ATTEMPT with the calculation of angles at every point, not so good

% 
% % erode slightly the cells so that small peaks along the boundary are eroded
% dataIn_L                = bwlabel(imopen(dataIn,strel('disk',1)));
% %
% 
% areaCells               = regionprops(dataIn_L,'Area');
% 
% 
% filledCells             = bwlabel(ismember(dataIn_L,find([areaCells.Area]>250)));
% 
% numCells                = max(filledCells(:));
% 
% stepBoundary            = 7;
% distanceFromEdge        = 6;
% anglePointyBit          = 115*pi/180;
% figure(1)
% hold off
% imagesc(filledCells)
% hold on
% clear finalPointBits
% %%
% for counterCells =1:numCells
%     %counterCells =1;
%     finalPointBits{counterCells} =[];
%     finalPoints{counterCells} =[];
%     %%
%     %stepBoundary            = 7;
%     [initialR,initialC] = find(filledCells==counterCells);
%         % % Calculate the radius from a roughly calculated centre
%     meanCurrentCell         = mean([initialR,initialC]);
%     currentContour          =(bwboundaries(filledCells==counterCells));
%     rhoContour1             = ([currentContour{1}(:,2)-meanCurrentCell(2) currentContour{1}(:,1)-meanCurrentCell(1)]);
%     rhoContour2             = rhoContour1.*rhoContour1;
%     rhoContour              = sqrt(sum(rhoContour2,2));
%     [q1,q2]                 =min(rhoContour);
%     initialPoint            =(currentContour{1}(q2,:));
%     %initialPoint =  [157 57];
%     %initialPoint =  [28 108];
%     %initialPoint =  [71 151];
%     
%     %%
%     % Find all the points that define the boundary of the cell
%     contourCell         = bwtraceboundary(filledCells,initialPoint, 'n', 8,inf,'clockwise');
%     numPointsBoundary   = size(contourCell,1);
%     for stepBoundary=5:9
%         % Subsample the points so that the slopes are more significant and noise is decreased
%         % it is important to append the first points at the end to guarantee circular detection
%         %contourCell1        = [contourCell(5:stepBoundary:end,:);contourCell([1 stepBoundary],:)];
%         contourCell1        = [contourCell(5:stepBoundary:end,:)];
%         
%         
%         % Calculate angles between the sampled points of the boundary
%         contourCell2        = diff(contourCell1(:,2))+1i*diff(contourCell1(:,1));
%         angleContour        = (angle(contourCell2 ));
%         
%         % unwrap the angle to avoid having differences between points that are close to -pi,+pi
%         % which are in reality very close to each other
%         angleContourUnwrap  = unwrap(angleContour);
%         
%         %find the difference between neighbouring angles and if they exceed a predefined value
%         %(115 degrees approx) determine it is a pointy edge
%         angleContourDiff    = diff(angleContourUnwrap);
%         numPointyBits       = find(angleContourDiff>anglePointyBit);
%         
% %             figure     
% %         hold off
% %            plot((angleContour),'b-','LineWidth',1);
% %         hold on
% %            plot(angleContourUnwrap,'r','LineWidth',1);
% %            plot(angleContourDiff,'m-o','LineWidth',2);
% %     
% %            grid on
%         
%         % angleContour2       = 0.5*angleContour(1:end-1)+0.5*angleContour(2:end);
%         % angleContour3       = 0.33*angleContour(1:end-2)+0.33*angleContour(2:end-1)+0.33*angleContour(3:end);
%         % angleContour4       = 0.25*angleContour(1:end-3)+0.25*angleContour(2:end-2)+0.25*angleContour(3:end-1)+0.25*angleContour(4:end);
%         
%         %
%         %angleContourFilt    = imfilter(angleContour,gaussF(1,1,1)','replicate');
%         
%         %plot(contourCell(:,2),contourCell(:,1),'g','LineWidth',2);
%         if ~isempty(numPointyBits)
%             numPointyBits2       =rem(numPointyBits*stepBoundary,numPointsBoundary);
%             numPointyBits2(numPointyBits2==0)=1;
%             disp([counterCells stepBoundary numPointyBits2'])
%             % if the points are far from the previously detected ones, add to the register
%             %checkDistance                    =(abs(finalPointBits{counterCells}-numPointyBits2))<9;
%             %if ~all(checkDistance)
%             finalPoints{counterCells}           = [finalPoints{counterCells}; contourCell1(numPointyBits,:)];
%             finalPointBits{counterCells}        = [finalPointBits{counterCells} numPointyBits2'];
%             %end
%             
%             
%             plot(contourCell1(numPointyBits,2),contourCell1(numPointyBits,1),'wx','LineWidth',2);
%         end
%         
%     end
%     if ~isempty (finalPointBits{counterCells})
%         %Finally remove any points in the boundary
%        closeTopLeft = any(finalPoints{counterCells}<=distanceFromEdge,2);
%        closeRight   = any(finalPoints{counterCells}(:,1)>=rows-distanceFromEdge+1,2);
%        closeBottom  = any(finalPoints{counterCells}(:,2)>=cols-distanceFromEdge+1,2);
%        finalPointBits{counterCells}(closeTopLeft|closeRight|closeBottom)=[];
%        %finalPointBits{counterCells}()=[];
%        %finalPointBits{counterCells}()=[];
%        
%        finalPoints{counterCells}(any(finalPoints{counterCells}<=5,2),:)=[];
%        finalPoints{counterCells}(any(finalPoints{counterCells}(:,1)>=rows-4,2),:)=[];
%        finalPoints{counterCells}(any(finalPoints{counterCells}(:,2)>=cols-4,2),:)=[];
%        
%        % remove duplicate points that are very close to one another
%         %finalPointBits{counterCells}        = sort(finalPointBits{counterCells});
%         if ~isempty(finalPointBits{counterCells})
%             [q1,q2]                             = sort(finalPointBits{counterCells});
%             PointsToKeep                        =q2([1 find(diff((q1))>9)+1]);
%             
%             finalPointBits{counterCells}        = finalPointBits{counterCells}(PointsToKeep);
%             finalPoints{counterCells}           = finalPoints{counterCells}(PointsToKeep,:);
%             % check circular condition
%             
%             if (finalPointBits{counterCells}(1)+numPointsBoundary - finalPointBits{counterCells}(end))<9
%                 finalPointBits{counterCells}= finalPointBits{counterCells}(1:end-1);
%             end
%             
%             
%             plot(finalPoints{counterCells}(:,2),finalPoints{counterCells}(:,1),'ro','LineWidth',2);
%         end
%        %plot(contourCell1(numPointyBits,2),contourCell1(numPointyBits,1),'ro','LineWidth',2);
%     end
%     %%
% %         figure(2)
% %         hold off
% %            plot((angleContour),'b-','LineWidth',1);
% %         hold on
% %            plot(angleContourUnwrap,'r','LineWidth',1);
% %            plot(angleContourDiff,'m-o','LineWidth',2);
% %     
% %            grid on
% %          figure(4)
% %     
% %          polar(angleContour,rhoContour(2:end),'k-o')
% %%
% end













%%




% 
% % Obtain Skeletons, when pointy, a branch of the skeleton should traverse the point
% dataInSkel              = bwmorph(dataIn,'skel',inf);
% dataInThin              = bwmorph(dataIn,'thin',inf);
% sizeToErode             = 7;
% % erode the data so that only the edges of the skeleton remain
% %dataInEroded            = imopen(imerode(dataIn,strel('disk',5)),ones(1));
% dataInEroded            = (imerode(dataIn,strel('disk',6)));
% pointyBits              = dataInSkel.*(dataIn-dataInEroded);
% dataInSkel2              = bwmorph(dataInEroded,'skel',inf);
% dataInThin2              = bwmorph(dataInEroded,'thin',inf);
% dataInThin3              = bwmorph(dataIn,'thin',sizeToErode);
% %
% figure(1)
% imagesc((1-dataInEroded).*dataInThin3)
% figure(2)
% imagesc((dataInEroded)+2*dataInThin3)
% 
% 
% %%
% sizeThin =7;
% dataInOpened            = imerode(dataIn,ones(sizeThin));
% dataInEroded            = imerode(dataIn,ones(sizeThin));
% dataInEroded3            = imopen(dataIn,ones(3));
% dataInThinned           = bwmorph(dataIn,'thin',sizeThin);
% imagesc(+dataInEroded3-dataInEroded)
% 
% 
% %%
% sizeToErode             = 7;
% dataInThin3              = bwmorph(dataIn,'thin',sizeToErode);
% dataInEroded            = imopen(imerode(dataIn,strel('disk',sizeToErode-2)),ones(1));
% %%
% VlargeSize = 19;
% largeSize  = 10;
% allPointyBits_L         = bwlabel((1-dataInEroded).*dataInThin3);
% allPointyBits_P         = regionprops(allPointyBits_L,'Area');       
% VlongPointyB             = ismember(allPointyBits_L,find([allPointyBits_P.Area]>VlargeSize));
% longPointyB             = ismember(allPointyBits_L,find(([allPointyBits_P.Area]>largeSize).*([allPointyBits_P.Area]<=VlargeSize)));
% smallPointyB             = ismember(allPointyBits_L,find([allPointyBits_P.Area]<=largeSize));
% 
% %figure(3)
% imagesc((dataIn)+4*smallPointyB+5*longPointyB+6*VlongPointyB)
% 
% %%
% % of the points that remain, label to find the axes lengths, 
% % look a) for very long and b) long and thin
% pointyBits_L            = bwlabel(pointyBits);
% pointyBits_A            = regionprops(pointyBits_L,'MajorAxisLength','MinorAxisLength');
% largePointyB            = ismember(pointyBits_L,find([pointyBits_A.MajorAxisLength]>20));
% largeThinPointyB        = ismember(pointyBits_L,find((([pointyBits_A.MajorAxisLength]<=20).*[pointyBits_A.MajorAxisLength]>12).*([pointyBits_A.MinorAxisLength]<5)));
% 
% finalPointyBits         = largePointyB+largeThinPointyB;
% %imagesc(dataIn+dataInSkel+1*largePointyB+1*largeThinPointyB)
