function [G] = RPBComputeGradient(k, grad_raster_time, gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %Constant
    G = zeros(size(k));
    for idx=2:size(k,2) %number points in trajectory
        G(:,idx) = (k(:,idx) - k(:,idx-1))/ (gamma*grad_raster_time);
    end
end

