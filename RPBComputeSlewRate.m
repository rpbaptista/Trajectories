function [slew_rate] = RPBComputeSlewRate(p,dx,nbPoints)

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dx','var') || isempty(dx)
  dx = 5e-3;
end
if ~exist('nbPoints','var') || isempty(nbPoints)
  nbPoints = 672;
end
% Constants
%nbPoints = 672;
dt = 10e-6;
gamma = 11.26e6;
T2e = nbPoints * dt;
L = 180e-3;
Kmax = 1/(2*dx);

G = 1/(2*gamma*T2e*dx) * (1+2*p.^3)./(3*p.^2) ;
slew_rate = (gamma*G.^2*L)./p;
%slew_rate = (1+2*p.^3)./(3*p.^2)
end

%L = 1/(Kmax*sin(atan((1/FOV)/(1/Kmax))))
