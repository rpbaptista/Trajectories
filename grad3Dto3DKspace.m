function [kspace] = grad3Dto3DKspace(grad3D, grad_raster_time, gamma)
        kx = RPBComputeTrajectory(squeeze(grad3D(1,:,:))', grad_raster_time, gamma);
        ky = RPBComputeTrajectory(squeeze(grad3D(2,:,:))', grad_raster_time, gamma);
        kz = RPBComputeTrajectory(squeeze(grad3D(3,:,:))', grad_raster_time, gamma);
        kspace = permute(cat(3, kx, ky, kz),[3,2,1]);

        
end

