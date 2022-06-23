function [] = RPBWriteGradientFile(outputfile,...
                                                        kSpace, grad3D)
    
    %RPBwriteGradientFile This functions writes .bin compatible with CSGRE
    % sequence based on 
   
    %% Verifying input arguments
    [~,~,ext] = fileparts(outputfile);
    if ~strcmp(ext,'.bin'); disp('Outputfile should be .bin'); return; end     

    % BIN file to write 
    fileoutID = fopen(outputfile,'w');
    if exist('ResultFile.txt', 'file')==2; delete(outputfile); end
    if fileoutID == -1; disp('Couldnt open for write the file'); return; end
    
    %% Writing .bin
    % Getting dimensions
    [dimensions, NumOfSamples, NumOfSegments] = size(kSpace);
    fprintf('Dimensions: %d\nNumOfSamples: %d\nNumOfSegments: %d\n',dimensions,NumOfSamples, NumOfSegments);
    
    % Dimension; number of segments; num of samples
    fwrite(fileoutID,dimensions,'float');
    fwrite(fileoutID,NumOfSegments,'float');
    fwrite(fileoutID,NumOfSamples,'float');
  
    % Writing KStarts
    kStart = kSpace(:,1,:);
    kStart = kStart(:);
    fwrite(fileoutID,kStart,'float');

    % Reading gradient values to find the maximum for each direction
    gMax = max(abs(grad3D(:)));
    
    % Writing max
    fwrite(fileoutID, single(1000*gMax),'float');
    fprintf('GradientMax: %.02d mT\n',1000*gMax);
    
    % Normalizing and writing
    grad3D = grad3D/gMax;
    for idx_projection = 1:NumOfSegments
        grad_project = grad3D(:,:,idx_projection) ;
        fwrite(fileoutID, grad_project(:),'float');
    end

    %% Close file
    fclose(fileoutID);

end

