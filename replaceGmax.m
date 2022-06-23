function  replaceGmax(filename_bin, ratio)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 %% Verifying input arguments
 
    [grad, gradMax, kStarts ] = readGradientBinFile(filename_bin);
    max(ratio*gradMax*grad(:))
   % gradMax = 1;
   grad_out = ratio*gradMax*grad/1e3;
   [flag_viability_x, ~] = RPBCalcCheckSlewRate(grad_out(1,:,:));
   [flag_viability_y, ~] = RPBCalcCheckSlewRate(grad_out(2,:,:));
   [flag_viability_z, ~] = RPBCalcCheckSlewRate(grad_out(3,:,:)); 
   if (flag_viability_x && flag_viability_y && flag_viability_z)
     RPBWriteGradientFile('fake_2.bin',zeros( 3,1055,10404),grad_out);
   else
       disp('shit')
   end
end

