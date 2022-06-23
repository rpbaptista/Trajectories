% read trajetory comùpute psf

path_bin = {
    'Y:\Renata_Baptista\Trajectories\2020-09-21\dar_1056s_18000s_3mm_100p_23Na_RPB.bin';
    'Y:\Renata_Baptista\Trajectories\2020-09-21\tpi_1056s_18000s_3mm_40p_23Na_RPB.bin';};
names = {
    'radial.mat',
    'tpi.mat'}
for n = 1:numel(path_bin)
    trajectory = readGradientBinFile(path_bin{n});
    trajectory = permute(trajectory(:,:,:),[2,3,1]);
    save(names{n},'trajectory');
  %  kx = permute(trajectory(1,:,:),[3,2,1]);
  %  ky =  permute(trajectory(2,:,:),[3,2,1]);
  %  kz =  permute(trajectory(3,:,:),[3,2,1]);
  %  ksp = [kx;ky;kz];
  %  fft_kspace = psf(trajectory)
end