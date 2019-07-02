function [dataIn]=readMTIFF_G_R(dataInName)
% function [dataIn]=readMTIFF(dataInName)
%
% Read a multi-frame TIFF file and save as a single Matlab matrix called
% dataIn
%
% It is assumed that the channels alternate one GREEN one Red

%% Parse input
switch nargin
    case 0
        %----- no data received, Open question dialog and pass to next section to analyse
        [pathname]                          =  uigetfile('*.tiff','Please select Multiple Tiff File');
        if pathname~=  0
            % pass the pathname to same function to process
            [dataIn]                        = readMTIFF_G_R(pathname);
        else
            %disp('Folder not found');
            dataIn=[];
            return;
        end
    case 1
        %----- one argument received, it can be a char of a matlab name or a folder
        if isa(dataInName,'char')
            % dataInName is a file name, can be .mat or a directory
            if (strcmp(dataInName(end-2:end),'tif'))||(strcmp(dataInName(end-3:end),'tiff'))
                try
                    numImages                   = numel(imfinfo(dataInName));
                    % Read first for dimensions
                    firstImage                  = imread(dataInName,1);
                    [rows,cols,levs]            = size(firstImage);
                    % initialise dataIn
                    dataIn(rows,cols,3,numImages/2) = uint8(0);
                    for counterImages =1:2:numImages-1
                        disp(strcat('Reading image',32,num2str(counterImages),32,'/',32,num2str(numImages)))
                        positionTime = floor((counterImages+1)/2);
                        dataIn(:,:,2,positionTime) = imread(dataInName,counterImages);
                        dataIn(:,:,1,positionTime) = imread(dataInName,counterImages+1);
                        
                    end
                catch
                    disp('Could not read Multiple TIFF file')
                end
            else
                disp('Could not read Multiple TIFF file')
                dataIn=[];
            end
        else
            dataIn=[];
        end
end
end
