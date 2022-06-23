FOV = 180e-3;
dx = 5e-3;
theta = atan((1/FOV)/(1/2*dx));
[kx, ky, kz] = RPBgeneratePointsFromTrajectory(theta,0);

gx = RPBComputeGradient(kx);
gy = RPBComputeGradient(ky);
gz = RPBComputeGradient(kz);
[valididity_tag_x, sr_x]= RPBCalcCheckSlewRate(gx);
[valididity_tag_y, sr_y]= RPBCalcCheckSlewRate(gy);
[valididity_tag_z, sr_z]= RPBCalcCheckSlewRate(gz);
%%
figure(1)
plot3(kx,ky,kz)

[kx, ky, kz] = RPBgeneratePointsFromTrajectory(30,30);
figure(2)
plot3(kx,ky,kz)


