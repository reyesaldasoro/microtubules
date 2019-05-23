function [dataIn]=readMTIFF(dataInName)
% function [dataIn]=readMTIFF(dataInName)
%
% Read a multi-frame TIFF file and save as a single Matlab matrix called
% dataIn

%% Parse input
switch nargin
    case 0
        %----- no data received, Open question dialog and pass to next section to analyse
        [pathname]                          =  uigetfile('*.tiff','Please select Multiple Tiff File');
        if pathname~=  0
            % pass the pathname to same function to process
            [dataIn]                        = readMTIFF(pathname);
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
                    for counterImages =1:numImages
                        disp(strcat('Reading image',32,num2str(counterImages),32,'/',32,num2str(numImages)))
                        dataIn(:,:,:,counterImages) = imread(dataInName,counterImages);
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
