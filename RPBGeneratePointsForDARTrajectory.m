% Generate Points for Radial Density Adapted Trajectory
function [kx,ky,kz] = RPBGeneratePointsForDARTrajectory(theta,phi,...
                                                      p, spatialResolution,...
                                                      nPointsKspaceTrajectory,...
                                                      gamma, dt, grad_raster_time)
                                                  
    %GeneratePointsForTrajectory It computes kx,ky,kz for all angles for
    %3DRADial Density
    %sequence
    % p = 1 -> radial

    %% Parameters
    % CONSTANTS
    dwellTime = 10e-6;

    
    %  Treating missing arguments
    if ~exist('gamma','var');   gamma = 11.26e6; end; %sodium
    if ~exist('limit_slew_rate','var');   limit_slew_rate = 180-3; end; %sodium
    if ~exist('dt','var');      dt = 10e-6;      end;
    if ~exist('grad_raster_time','var') ; grad_raster_time = 10e-6; end;
    fprintf('-----Parameters-------:\np %d\nGamma: %d\ndt: %d\nGradRaster%d\n\n',p,gamma, dt, grad_raster_time);
    
   %% Computing acessory variables 

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    % avoid numerical errors
    cos_theta(abs(cos(theta))< 1.0e-7) = 0;
    
    % Computing values depending on the inputs 
    kmax = 1 /(2*spatialResolution); 
    
    requiredGradient = (1/(2*gamma*nPointsKspaceTrajectory*dt*spatialResolution))*((1+2*p^3)/(3*p^2)); 
    %requiredGradient = (1/(2*gamma*nPointsKspaceTrajectory*dt*spatialResolution)); 
    %k0 = 2*gamma*requiredGradient^2.0/limit_slew_rate;
    k0 = p*kmax;
    
    rampTime = 20*grad_raster_time;
    rampSteps = (floor(rampTime/dwellTime));
    rampMoment = 0.5*gamma*requiredGradient*double(rampSteps)*dwellTime;
    plateauSteps = floor(((k0 - double(rampMoment))/(gamma*requiredGradient))/dwellTime);
    if plateauSteps + rampSteps > nPointsKspaceTrajectory
        plateauSteps = nPointsKspaceTrajectory - rampSteps;
    end
    RadialPoints = max(nPointsKspaceTrajectory - double(plateauSteps) - double(rampSteps),0 ); % MORE POINTS TO COMPUTE GRADIENT
    
    %   Initializing outputs
    kx = zeros(size(theta,1),nPointsKspaceTrajectory);
    ky = zeros(size(theta,1),nPointsKspaceTrajectory);
    kz = zeros(size(theta,1),nPointsKspaceTrajectory);
    %kx = zeros(size(theta,1),nPointsKspaceTrajectory);

    %% Computing the trajectory
    % ramp up
    i = 1;
    t = 0;
    for idx=1:rampSteps
        kx(:,i) = sin_theta*0.5*gamma*requiredGradient*t*t/(rampSteps*dwellTime);
        kz(:,i) = cos_theta*0.5*gamma*requiredGradient*t*t/(rampSteps*dwellTime);
        t = t + dwellTime;
        i = i + 1;
    end
    
    % plateau
    t = 0;
    if RadialPoints ~= 0
        for idx=1:plateauSteps+1% % Is it correct why 0 t + 1
            kx(:,i) =  sin_theta*(rampMoment + gamma*requiredGradient*t);
            kz(:,i) =  cos_theta*(rampMoment + gamma*requiredGradient*t);
            i = i + 1;
            t = t + dwellTime;
        end
    else
        for idx=1:plateauSteps% % Is it correct why 0 t + 1
            kx(:,i) =  sin_theta*(rampMoment + gamma*requiredGradient*t);
            kz(:,i) =  cos_theta*(rampMoment + gamma*requiredGradient*t);
            i = i + 1;
            t = t + dwellTime;
        end
    end
    
    % Radial t > t0
    actual_k0 = double(rampMoment + gamma*requiredGradient*plateauSteps*dwellTime);
    t = dwellTime;
   
    i_save = i;
   
    for idx=1:RadialPoints-1
        K_dar = (3*gamma*requiredGradient*(actual_k0^2.0)*t + actual_k0^3.0)^(1/3);
        
        kx(:,i) = K_dar.*sin_theta;  
        kz(:,i) = K_dar.*cos_theta;
        t = t + dwellTime;
        i = i + 1 ;
    end
   % Now change only in idx
    idx_special_case = find((theta == 0.0) | (abs(theta-pi)<1e-6));
    t = dwellTime;
    i = i_save;
    if RadialPoints ~= 0
        for idx=1:RadialPoints-1
            K_dar = (3*gamma*requiredGradient*(actual_k0^2.0)*t + actual_k0^3.0)^(1/3);
            kx(idx_special_case,i) = 0;
            ky(idx_special_case,i) = 0;
            kz(idx_special_case,i) = K_dar;
            t = t + dwellTime;
            i = i + 1 ;
        end
         
    end
    %rotate about z_axis
    for idx=1:nPointsKspaceTrajectory
        x = kx(:,idx);
        y = ky(:,idx);
        z = kz(:,idx);

        kx(:,idx) = x.*cos(phi) - y.*sin(phi);
        ky(:,idx) = x.*sin(phi) + y.*cos(phi);
        kz(:,idx) = z;
    end
end