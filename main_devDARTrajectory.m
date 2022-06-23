% Dev - Script to develop DAR trajectory
% Author : Renata Porciuncula Baptista
% e-mail : renata.porciunculabaptista@cea.fr


%%
[grad, gradMax, ~ ] = readGradientBinFile('Y:\Renata_Baptista\Trajectories\1H_sparkling_NS_102_cutoff_30_decay_3_Npts_1056_OS_2_FOV_240.bin');
grad = grad*gradMax;

[gradTPI, gradMaxTPI, ~ ] = readGradientBinFile('Y:\Renata_Baptista\Trajectories\tpi_1056s_21000s_3mm_40p_1H_RPB.bin');
gradTPI = gradTPI * gradMaxTPI;

%%
figure(1)
plot(sqrt(reshape(grad(1,:,1:20),1,[]).^2 + reshape(grad(2,:,1:20),1,[]).^2 + reshape(grad(3,:,1:20),1,[]).^2), '-g')
hold on
plot(sqrt(reshape(gradTPI(1,:,1:20),1,[]).^2 + reshape(gradTPI(2,:,1:20),1,[]).^2 + reshape(gradTPI(3,:,1:20),1,[]).^2), '-b')

%figure(2)
%plot(reshape(grad(2,:,1:20),1,[]), '-g')
%plot(reshape(gradSod(2,:,1:20),1,[]),'-b')

%figure(3)
%plot(reshape(grad(3,:,1:20),1,[]))
%plot(reshape(gradSod(3,:,1:20),1,[]))

%figure(4)
% plot3(kx(:,:),ky(:,:),kz(:,:),'o')
% for i=1:10%size(kx,1)
%     for j=1:size(kx,2)
%         plot3(kx(i,1:j),ky(i,1:j),kz(i,1:j),'o')
%         pause(0.001)
%     end
% end