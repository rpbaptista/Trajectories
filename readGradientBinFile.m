function [grad, gradMax, kStarts ] = readGradientBinFile(inputfile)

    %readGradientBinFile It reads gradients from a given .bin file

    %   Author: Renata Porciuncula Baptista
    %   E-mail: renata.porciunculabaptista@cea.fr

    if (nargin < 1)
        info = 'Please select binary file to read';
        [file,pathname]=uigetfile('*.bin',info);
        inputfile = fullfile(pathname, file);
    end
    
   fileID = fopen(inputfile,'r');
   
   dimensions = fread(fileID,1,'float');
   NumOfSegments = fread(fileID,1,'float');
   NumOfSamples = fread(fileID,1,'float');
   
   fprintf('Dimensions: %d\nNumOfSamples: %d\nNumOfSegments: %d\n',dimensions,NumOfSamples, NumOfSegments);

   kStarts = reshape(fread(fileID, dimensions*NumOfSegments,'float'), dimensions,NumOfSegments);
   
   gradMax = fread(fileID,1,'float');
   fprintf('GradientMax: %.02d\n',gradMax);

   grad = reshape(fread(fileID,dimensions*NumOfSamples*NumOfSegments,'float'),dimensions,NumOfSamples,NumOfSegments);
   
   fclose(fileID);
end

