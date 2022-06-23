% Author: Renata Porciuncula Baptista
% E-mail: re.porci.rb@gmail.com
% Date: 27 February 2020

% This script generates the trajectory based in a given set of parameters
% verify it's feability and writes the correspondent gradient file compa-
% tible to the sequence GRE arb

%% Parameters
% CONSTANTS

% VARIABLES
p = 0.3; % radial 3600 / radial, spiral
number_of_projections = 144;
spatialResolution = 4e-3;
nPointsKspaceTrajectory = 1248;
nuclei = '23Na';
trajectory = 'tpi'; %dar, rad also available

if strcmp(nuclei, '31P');  gamma = 17.23e6; end;% Hz/T 
if strcmp(nuclei, '23Na'); gamma = 11.26e6; end;% Hz/T 
if strcmp(nuclei, '1H');   gamma = 42.58e6; end;% Hz/T 

dt = 10e-6; % used in approximation of T2ef as Npt*dt
grad_raster_time = 10e-6; % conferi 
limit_slew_rate = 200e-3;
output_file = sprintf('./output/%s_%ds_%ds_%dmm_%dp_%s_RPB.bin',trajectory, ...
                                                         nPointsKspaceTrajectory,...
                                                         number_of_projections,...
                                                         1000*spatialResolution,...
                                                         p*100, nuclei);
%% Computing trajectory
[theta, phi] = RPBComputeAngleSandroVersion(number_of_projections);


[kx, ky, kz] = RPBGeneratePointsForXTrajectory(trajectory, theta, phi,...
                                    p, spatialResolution,...
                                    nPointsKspaceTrajectory, ...
                                    gamma, dt, grad_raster_time); 
    
kSpace = permute(cat(3, kx, ky, kz),[3,2,1]);

%%
gx = RPBComputeGradient(kx, grad_raster_time, gamma);
gy = RPBComputeGradient(ky, grad_raster_time, gamma);
gz = RPBComputeGradient(kz, grad_raster_time, gamma);

g3D = permute(cat(3, gx, gy, gz),[3,2,1]);
%% Verifying viability
[flag_viability_x, ~] = RPBCalcCheckSlewRate(gx,limit_slew_rate,dt);
[flag_viability_y, ~] = RPBCalcCheckSlewRate(gy,limit_slew_rate,dt);
[flag_viability_z, ~] = RPBCalcCheckSlewRate(gz,limit_slew_rate,dt);

%% Writing .bin
if (flag_viability_x && flag_viability_y && flag_viability_z)
    RPBWriteGradientFile(output_file,kSpace, g3D); %...                            
else
    disp('Error: Slew rates superior to the limit')
end

%% Debug reasons
save('bin.mat', 'g3D', 'kSpace', 'grad_raster_time', 'gamma')