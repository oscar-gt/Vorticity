close all; clear all; clc;

%  Oscar Garcia-Telles
%  AMATH 481 Hw 4
%  24 November 12

% ------------------------------ Movies --------------------------------
%
%  Here we use finite difference methods to solve
%  for vorticity and create 'movies' that show
%  its evolution over time. 
%

%  A, B, C matrices
%  A matrix
L = 10;    %  Domain bound
n = 64;    %  So matrix A will be dimensions (64^2 x 64^2)
xspan = (linspace(-L,L, n+1))';  %  Column vector of time span
dt = xspan(2) - xspan(1);        %  dt value. Using dx elsewhere

A = zeros(n^2, n^2);   %  Initializing matrix A
d = -4*ones(n,1);      %  Col vector of -4
D = diag(d);           %  Diagonal matrix with -4 on main diag
e1 = ones(n,1);

yy = sparse(n,n);
for i = [-(n-1) -1 1 (n-1)]         %  Building diagonals of 1's on an (64 x 64) block
    y1 = spdiags(e1, i, n,n);
    yy = yy + y1;
end

b1 = yy + D;  % (8 x 8) block of -4 on main diag and 1's on other diags
%  Adding blocks b1 to the matrix
for i = 1:n:n^2 - (n-1)
    A(i:i+(n-1),i:i+(n-1)) = b1;
end

%  Including diagonals of 1's
e2 = ones(n^2,1);
ii= spdiags([e2 e2], [-n n], n^2, n^2);
%spy(ii)
A = A + ii;
%  Including the I blocks on the upper right, lower left
A(1:n, n*n - (n-1):n*n) = eye(n,n);
A(n*n - (n-1):n*n, 1:n) = eye(n,n);

%  Changing A(1,1) to 2
A(1,1) = 2;
%A = full(A./(dt^2));                           %  Complete A matrix

%  Matrix B
b = zeros(n^2, n^2);
B = spdiags([e2 -e2 e2 -e2], [-((n^2)-n) -n n ((n^2)-n)], b);
B = B./(2*dt);                                 %  Complete B matrix
%B = full(B);

%  Matrix C
%  Constructing (8 x 8) block to build larger matrix
C = zeros(n^2, n^2);
cblock = spdiags([e1 -e1 e1 -e1], [-(n-1) -1 1 n-1], n, n);
%  Now constructing (64 x 64) C matrix using smaller block
for i = 1:n:(n*n)-n+1
    C(i:i+(n-1), i:i+(n-1)) = cblock;
end
C = C./(2*dt);
%C = full(C);                                  %  Complete C matrix

A = sparse(A);
B = sparse(B);
C = sparse(C);

%  Initializing 
Lx = 20;    % x and y domain
Ly = 20;
nx = n;     %  n should be 64 for this assignment
ny = n;
N = nx*ny;  %  N = 4096 for this assignment  
v = 0.001;

x2 = linspace(-Lx/2, Lx/2, nx + 1);  %  Account for periodicity
x = x2(1:nx);  %  x vector
y2 = linspace(-Ly/2, Ly/2, ny + 1);  %  Account for periodicity
y = y2(1:ny);  %  y vector

[X, Y] = meshgrid(x,y);           %  2D domain

%  Initializing FFT variables
% Spectral Values
kx = (2*pi/Lx) * [0:(nx/2 - 1) (-nx/2):-1];
ky = (2*pi/Ly) * [0:(ny/2 - 1) (-ny/2):-1];
%  Changing first values to avoid dividing by zero
kx(1) = 10^(-6); 
ky(1) = 10^(-6);

%  2D grid
[KX, KY] = meshgrid(kx, ky);
fft_factor = ((KX.^2) + (KY.^2));


% ----------------------- Making Different Movies ---------------------
%  Need:  
%  1. Two oppositely charged vortices next to each other
%  2. Two same charged vortices next to each other
%  3. Two pairs of oppositely charged vortices which can be made to collide
%     with each other
%  4. Random assortment, try 10-15 vortices.


%  Specifying different time span
tend = 12;
TT = 0:0.5:tend;  %  Longer time span
dx = TT(2) - TT(1);

%  Shifting parameters
a = 1;  %  Horizontal shift
b = 0;  %  Vertical Shift
%  Initial shape parameter
c = 1;
%  Amplitude
amp = 5;

%  Oppositely charged vortices
UL1 =  amp*exp((-(X - a).^2)./c - (((Y - b).^2)./20));
UL2 = -amp*exp((-(X + a).^2)./c - (((Y + b).^2)./20));
U1 = UL1 + UL2;
U1re = reshape(U1, N, 1);
%  Timestepping, solving for vorticity
[t,   W1] = ode45('funFFT',   TT, U1re, [], v, dx, A, B, C, fft_factor);

M1 = frames(W1, TT);
close all;
movie(M1,1,8);
close all;

%  2 Same charged vortices next to each other
u2 = initialw(X, Y, 1, -1, 0, 0, 1, 8, 8);
[t, W2] = ode45('funFFT', TT, u2, [], v, dx, A, B, C, fft_factor);
M2 = frames(W2, TT);
close all; 
movie(M2, 1, 8);
% 
% 
%  2 pairs of oppositely changed vortices colliding
%  Upward vortices
u3 =  initialw(X, Y, 1, -1, 8, 8, 1, 8, -8);
%  Downward vortices
u33 = initialw(X, Y, 1, -1, -8, -8, 1, -8, 8);
u03 = u3+u33;
[t, W3] = ode45('funFFT', TT, u03, [], v, dx, A,B,C, fft_factor);

 M3 = frames(W3, TT);
 close all;
 movie(M3, 1, 8);
 close all;
 
%  Random assortment of vortices

u41 =  initialw(X, Y, 1, -1, 3, 3, 1, 8, -8);
u42 = initialw(X, Y, 1, -1, -3, -3, 1, -8, 8);
u43 = initialw(X, Y, -6, -6, 1, -1, 1, 5, -5);
u44 = initialw(X, Y, 6, 6, 1, -1, 1, -8, 8);
u45 = initialw(X, Y, 6, -6, 6, -6, 1, 8, 8);
u46 = initialw(X, Y, -6, 6, 6, -6, 1, 8, 8);
u47 = initialw(X, Y, 3, -3, 3, -3, 1, 5, 5);
u48 = initialw(X, Y, -3, 3, 3, -3, 1, 5, 5);
u4 = u42+u41+u2+u43+u44+u45+u46+u47+u48;

[t, W4] = ode45('funFFT', TT, u4, [], v, dx, A, B, C, fft_factor);
M4 = frames(W4, TT);
close all;
movie(M4, 1, 8);
close all;

