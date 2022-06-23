function [kx,ky,kz] = RPBGeneratePointsForTPITrajectory(theta,phi,...
                                                      p, spatialResolution,...
                                                      nPointsKspaceTrajectory,...
                                                      gamma,dt, grad_raster_time)
                                                  
    %GeneratePointsForTrajectory It computes kx,ky,kz for all angles for TPI
    %sequence

    %% Parameters
    % CONSTANTS
    dwellTime = 10e-6;

    %  Treating missing arguments
    if ~exist('p','var');       p = 0.3;         end;
    if ~exist('gamma','var');   gamma = 11.26e6; end; %sodium
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
    initial_k0 = p*kmax;
  
    requiredGradient = (1/(2*gamma*nPointsKspaceTrajectory*dt*spatialResolution))*((1+2*p^3)/(3*p^2)); 
    rampTime = 20*grad_raster_time;
    rampSteps = (floor(rampTime/dwellTime));
    rampMoment = 0.5*gamma*requiredGradient*double(rampSteps)*dwellTime;
    plateauSteps = floor(((initial_k0 - double(rampMoment))/(gamma*requiredGradient))/dwellTime);
    if plateauSteps + rampSteps > nPointsKspaceTrajectory
        plateauSteps = nPointsKspaceTrajectory - rampSteps;
    end
    tpiPoints = max(nPointsKspaceTrajectory - double(plateauSteps) - double(rampSteps),0 ); % MORE POINTS TO COMPUTE GRADIENT
    
    %   Initializing outputs
    kx = zeros(size(theta,1),nPointsKspaceTrajectory);
    ky = zeros(size(theta,1),nPointsKspaceTrajectory);
    kz = zeros(size(theta,1),nPointsKspaceTrajectory);
    
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
    
     t = 0;
    if tpiPoints ~= 0
        for idx=1:plateauSteps+1% % Is it correct why 0 t + 1
            kx(:,i) =  sin_theta*(rampMoment + gamma*requiredGradient*t);
            kz(:,i) =  cos_theta*(rampMoment + gamma*requiredGradient*t);
            i = i + 1;
            t = t + dwellTime;
        end
    else
        for idx=1:plateauSteps+1% % Is it correct why 0 t + 1
            kx(:,i) =  sin_theta*(rampMoment + gamma*requiredGradient*t);
            kz(:,i) =  cos_theta*(rampMoment + gamma*requiredGradient*t);
            i = i + 1;
            t = t + dwellTime;
        end
    end
    
    % TPI
    K0 = (3*gamma*requiredGradient*initial_k0^2*dwellTime + initial_k0^3)^(1/3);
    chi0 = sqrt((K0^4/initial_k0^4)-1);
    phi0 = (chi0 + atan(1/chi0))./(2*sin_theta);
    
    initial_k0 = double(rampMoment + gamma*requiredGradient*plateauSteps*dwellTime);
    t = dwellTime;
   
    i_save = i;
   
    for idx=1:tpiPoints-1
        K_tpi = (3*gamma*requiredGradient*(initial_k0^2.0)*t + initial_k0^3.0)^(1/3);
        chi_tpi = sqrt(((K_tpi^4.0)/(initial_k0^4.0))-1);
        phi_tpi = ((chi_tpi + atan(1/chi_tpi))./(2*sin_theta))-phi0;

        kx(:,i) = K_tpi.*cos(phi_tpi).*sin_theta;
        ky(:,i) = K_tpi.*sin(phi_tpi).*sin_theta;
        kz(:,i) = K_tpi.*cos_theta;
        t = t + dwellTime;
        i = i + 1 ;
    end
   % Now change only in idx
    idx_special_case = find((theta == 0.0) | (abs(theta-pi)<1e-6));
    t = dwellTime;
    i = i_save;
    for idx=1:tpiPoints-1
        K_tpi = (3*gamma*requiredGradient*(initial_k0^2.0)*t + initial_k0^3.0)^(1/3);
        kx(idx_special_case,i) = 0;
        ky(idx_special_case,i) = 0;
        kz(idx_special_case,i) = K_tpi;
        t = t + dwellTime;
        i = i + 1 ;
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

