function [KspaceMag] = getMagnitude(Kspace)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 KspaceMag = squeeze(sqrt(Kspace(1,:,:).^2 +Kspace(2,:,:).^2+Kspace(3,:,:).^2 ));
 
end

