function [k] = RPBComputeTrajectory(G, grad_raster_time, gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    k = zeros(size(G));
    for idx=2:size(k,2) %number points in trajectory
         k(:,idx)= G(:,idx)*(gamma*grad_raster_time) + k(:,idx-1);
    end
end

