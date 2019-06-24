function  finalTubules2 = growLinearRegion(allTubules)
% Process all tubules,
%  find maximum and proceed
%     connect neighbour that presents the least change in intensity
%     proceed until all neighbours are zero
%  find next maximum and iterate
%%
% Only take into account regions whose linear intensity is above 1.5
currentIntensities                  = allTubules.*(allTubules>1.5);

%delete the edges to avoid boundary conditions
currentIntensities([1 end],:)       =0;
currentIntensities(:,[1 end])       =0;
%%

currentTubule                       = 0;
finalTubules                        = zeros(size(allTubules));

discardCentreNeigh                  = [0 0 0;0 inf 0 ;0 0 0];


% Start iteration over the tubules
%%
while (max(currentIntensities(:)>3))
%%    
    currentTubule                   = currentTubule+1;
    %disp(currentTubule)
    %beep
    %pause(2)
    % find maximum
    [maxRow,maxCol]                 = find(currentIntensities==max(currentIntensities(:)));
%    
    %[maxPerCol,indMaxCol]       = max(currentIntensities);
    %[absMax,maxCol]             = max(maxPerCol);
    %absMaxF                     = absMax;
    %maxRow                      = indMaxCol(maxCol);
    frontPixelR                     = (  maxRow(1));
    frontPixelC                     = ( maxCol(1) );
    currentFront                    = currentIntensities(frontPixelR,frontPixelC);
    frontNeigh                      = abs(currentFront -currentIntensities(frontPixelR-1:frontPixelR+1,frontPixelC-1:frontPixelC+1))+discardCentreNeigh;
    
    [newRow,newCol]                 = find(frontNeigh==min(frontNeigh(:)));
    endPixelR                       = (  frontPixelR+newRow(1)-2);
    endPixelC                       = ( frontPixelC+newCol(1)-2 );
    currentEnd                      = currentIntensities(endPixelR,endPixelC);
    
    
    
    finalTubules(frontPixelR,frontPixelC)       = currentTubule;
    finalTubules(endPixelR,endPixelC)           = currentTubule;
    currentIntensities(frontPixelR,frontPixelC) = -1;
    currentIntensities(endPixelR,endPixelC)     = -1;
    
    endNeigh                        = abs(currentFront -currentIntensities(endPixelR-1:endPixelR+1,endPixelC-1:endPixelC+1))+discardCentreNeigh;
    frontNeigh                      = abs(currentFront -currentIntensities(frontPixelR-1:frontPixelR+1,frontPixelC-1:frontPixelC+1))+discardCentreNeigh;
%    
    %twoNeighs =[frontNeigh endNeigh];
    
    % Start iteration to grow  the current tubule
    %while (any(frontNeigh(:)<currentFront))&(any(endNeigh(:)<currentEnd))
    while any([((frontNeigh(:)<currentFront)); ((endNeigh(:)<currentEnd))])
%
        minEnd                      = min(endNeigh(:));
        minFront                    = min(frontNeigh(:));
        
        %disp([minEnd minFront currentFront])
        
        if (minFront<=minEnd)
            % Pixel to join becomes the front
            [newRow,newCol]         = find(frontNeigh==min(frontNeigh(:)));
            
            % move the position of the front pixel
            prevPixR                = frontPixelR;
            prevPixC                = frontPixelC;
            
            frontPixelR             = ( frontPixelR+newRow(1)-2);
            frontPixelC             = ( frontPixelC+newCol(1)-2 );
            
            %remove previous pixel from the data
            finalTubules(frontPixelR,frontPixelC)       = currentTubule;
            currentIntensities(frontPixelR,frontPixelC) = -1;
            
        else
            % Pixel to join becomes the front
            [newRow,newCol]         = find(endNeigh==min(endNeigh(:)));
            
            % move the position of the front pixel
            prevPixR                = endPixelR;
            prevPixC                = endPixelC;
            
            endPixelR               = ( endPixelR+newRow(1)-2);
            endPixelC               = ( endPixelC+newCol(1)-2 );
            
            %remove previous pixel from the data
            finalTubules(endPixelR,endPixelC)       = currentTubule;
            currentIntensities(endPixelR,endPixelC) = -1;
        end
%        
        % discard corners and neighbours
        toDiscardV                  = currentIntensities(prevPixR-1:prevPixR+1,prevPixC+newCol(1)-2);
        toDiscardH                  = currentIntensities(prevPixR+newRow(1)-2,prevPixC-1:prevPixC+1);
        toDiscardV(toDiscardV>0)    = -2;
        toDiscardH(toDiscardH>0)    = -2;
        
        currentIntensities(prevPixR-1:prevPixR+1,prevPixC+newCol(1)-2) =toDiscardV;
        currentIntensities(prevPixR+newRow(1)-2,prevPixC-1:prevPixC+1) =toDiscardH;
        currentIntensities(prevPixR,prevPixC) = -1;
        %currentEnd                 = currentIntensities(endPixelR,endPixelC);
        currentEnd                  = mean(allTubules(currentIntensities==-1));
        currentFront                = currentEnd;
        
        %recalculate neighbourhoods
        endNeigh                    = abs(currentEnd -currentIntensities(endPixelR-1:endPixelR+1,endPixelC-1:endPixelC+1))+discardCentreNeigh;
        frontNeigh                  = abs(currentFront -currentIntensities(frontPixelR-1:frontPixelR+1,frontPixelC-1:frontPixelC+1))+discardCentreNeigh;
        %imagesc(currentIntensities)
        %drawnow;
    end
    % dilate the current tubule and discard all neighbours to avoid having any contact to
    % it
    neighboursToDelete = (currentIntensities.*(bwmorph(finalTubules,'dilate')-finalTubules))>0;
    currentIntensities(neighboursToDelete) = -2;
       % imagesc(finalTubules)
        %drawnow;
%%        
end

%% Clean the microtubules
finalTubules2                       = bwmorph(finalTubules,'spur',1);
finalTubules_R                      = regionprops(finalTubules,'Area','MajorAxisLength','MinorAxisLength');

shortTubules                        = ismember(finalTubules,find([finalTubules_R.Area]<9));
mediumTubules                       = ismember(finalTubules,find([finalTubules_R.Area]<20));
notElongatedTubules                 = ismember(finalTubules,find(([finalTubules_R.MajorAxisLength]./[finalTubules_R.MinorAxisLength])<2.75));

finalTubules2                       =finalTubules.*(1-((notElongatedTubules).*(mediumTubules)+(shortTubules))>0);

%imagesc(finalTubules2)

