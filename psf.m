function PSF = psf(ksp)


%z = pupil.*exp(2*pi*j*w);
%z = fftshift(z);
%z = fft2(z);
%z = fftshift(z);
%z = z.*conj(z);
%zmax = max(max(z));
%z = z/zmax;
    ksp = complex( zeros( szRead, szY, szZ, NumSamplings, 'single' ) );

    ksp = permute( ksp, [ 2 1 3 4 ] );

    % Calculate and normalize the Point Spread Functions
    %
%    shifts = [0  0  0  0 ] ./ 2;
    PSF = abs( ifft(ifft(ifft(  ksp, [],1),[],2),[],3) );
 %   PSF = cmshiftnd( PSF, shifts );
    PSF = bsxfun( @rdivide, PSF, max(max(max( PSF,[],1),[],2),[],3) );
