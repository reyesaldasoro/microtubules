function [finalTubules,res_strength,res_orientation]=linearStructures(dataIn,cellBody)

% find dimensions of the input data
[rows,cols]                                     = size(dataIn);

%define all the orientations that will be used to detect linear structures
orientation{1}                                  = [3  8  13 18 23];         % 0     degrees
orientation{2}                                  = [4  8  9  13 17 18 22];   % 22.5  degrees
orientation{3}                                  = [5  9  13 17 21];         % 45    degrees
orientation{4}                                  = [10 9  14 13 12 17 16];   % 77.5  degrees
orientation{5}                                  = [11 12 13 14 15];         % 90    degrees
orientation{6}                                  = [6  7  12 13 14 19 20];   % 112.5 degrees
orientation{7}                                  = [1  7  13 19 25];         % 135   degrees
orientation{8}                                  = [2  7  8  13 18 19 24];   % 157.5 degrees

%fieldOne       will contain the elements on one side of the orientation line
%fieldTwo       will contain the elements on the other side of the orientation line
%fieldOutside   will contain the elements on one side of the orientation line

fieldOne{1}                                     = [1  2  6  7  11 12 16 17 21 22];
fieldOne{2}                                     = [1  2  3  6  7  11 12 16 21   ];
fieldOne{3}                                     = [1  2  3  4  6  7  8  11 12 16];
fieldOne{4}                                     = [1  2  3  4  5  6  7  8  11   ];
fieldOne{5}                                     = [1  2  3  4  5  6  7  8  9  10];
fieldOne{6}                                     = [1  2  3  4  5  8  9  10 15   ];
fieldOne{7}                                     = [2  3  4  5  8  9  10 14 15 20];
fieldOne{8}                                     = [3  4  5  9  10 14 15 20 25   ];

fieldOutside{8}                                 =[];
fieldTwo{8}                                     =[];


elemsOrientation                                ={5, 7, 5, 7, 5, 7, 5, 7};
elemsFieldOne                                   ={10,9,10, 9,10, 9,10, 9}; 

for kk=1:8
    fieldOutside{kk}                            = setdiff((1:25),orientation{kk});
    fieldTwo{kk}                                = setdiff(fieldOutside{kk},fieldOne{kk});
end

nn=2;


res_strength(rows,cols)                         = 0;
res_orientation(rows,cols)                      = 0;


for counterR = 1+nn:1:rows-nn
    %disp(counterR)
    for counterC = 1+nn:1:cols-nn
        
        tempdata                                = double(dataIn(counterR-nn:counterR+nn,counterC-nn:counterC+nn));
        %clear lineStrength
        lineStrength                                    = zeros(1,8);
        for kk=1:8
            %kk=1;
            
            
            %tempDataOrientation1 = (tempdata(orientation{kk}));
            tempDataOrientation                 = sum(tempdata(orientation{kk}))/elemsOrientation{kk};
            %tempDataOrientation = mean(tempDataOrientation1);
%            tempDataOutside     = mean(tempdata(fieldOutside{kk}));

            %tempDataOne1        = (tempdata(fieldOne{kk}));
            tempDataOne                         = sum(tempdata(fieldOne{kk}))/elemsFieldOne{kk};            
            %tempDataOne         = mean(tempdata(fieldOne{kk}));
            %tempDataTwo1        = (tempdata(fieldTwo{kk}));
            tempDataTwo                         = sum(tempdata(fieldTwo{kk}))/elemsFieldOne{kk};            
            
            %tempDataTwo         = mean(tempdata(fieldTwo{kk}));
            
            %[mean(tempDataOrientation)-mean(tempDataOutside) mean(tempDataOne) mean(tempDataTwo) mean(tempDataOrientation)]
            
            tempDataOutside                     = (tempDataOne+tempDataTwo)/2;
            %tempInter1                          = (tempDataOrientation-tempDataOutside);
            %tempInter2                          = (tempDataOrientation>tempDataOne);
            %tempInter3                          = (tempDataOrientation>tempDataTwo);
            %lineStrength (kk)                   = tempInter1.*(tempInter2&tempInter3);
            
            lineStrength (kk)                   = (tempDataOrientation-tempDataOutside).*(tempDataOrientation>tempDataOne).*(tempDataOrientation>tempDataTwo);
        end
        %tempdata2=tempdata;
        %     tempdata2(orientation{kk})=-1;
        %
        %imagesc(tempdata2)
        %drawnow
        %if any(lineStrength>0)
        [finalLineStrength,finalOrientation]    = max(lineStrength);
        
        res_strength(counterR,counterC)         = finalLineStrength;
        res_orientation(counterR,counterC)      = finalOrientation;
        
    end
end

if exist('cellBody','var')
%%    
    cellBodyD                                   = imdilate(cellBody,strel('disk',5));
    allTubules                                  = res_strength.*(1-cellBodyD);
   % strongTubules       = ((res_strength>3.5).*(1-cellBodyD));
   % weakTubules         = imclose(bwlabel((res_strength>2.5).*(1-cellBodyD)),ones(2));
   %  weakNextStrong      = unique(strongTubules.*(weakTubules));
    finalTubules                                = growLinearRegion(allTubules);
    % Process tubules as linear structures in which each element can only have 2
    % neighbours, region-grow from maximum until all new neighbours are below a level
    
    
%     % keep weak tubules with a peak over the strong ones
%     figure(22)
%     hold on
%     [yyy,xxx] = hist(res_strength(:).*(1-cellBodyD(:)),100);
%     plot(xxx,yyy)
%     finalTubules_1      = ismember(weakTubules,weakNextStrong(2:end));
%     finalTubules_2      = bwmorph(finalTubules_1,'thin','inf');
%     finalTubules_3      = bwmorph(finalTubules_2,'spur',1);
%     finalTubules        = bwmorph(finalTubules_3,'clean');
    
   %imagesc(double(dataIn)+50*(finalTubules>0));
%%    
else
    finalTubules                                = zeros(rows,cols);
    
end

%[finalLineStrength finalOrientation]