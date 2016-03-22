function r = funFFT(tspan,u,dummy,k,dx,A, B, C, fft_factor)

%  Changing vector into a matrix
umat = reshape(u, 64,64);
%  FFT
ut = fft2(umat);
utk = -ut ./ fft_factor;
utkinv = real(ifft2(utk));
psi = reshape(utkinv,4096, 1);


%  Time stepping
r = (k.*A*u) - (B*psi).*(C*u) + (C*psi).*(B*u);

end