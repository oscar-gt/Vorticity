%  Calculates the initial vorticity with spacial domain X,Y, and spacial
%  distortion parameters a,b,c. Amp1 is the amplitude of one gaussian bump,
%  amp2 to for the second
function u = initialw(X, Y, a, aa, b, bb, c, amp1, amp2)

UL1 = amp1*exp((-(X + a).^2)./c - (((Y + b).^2)./20));
UL2 = amp2*exp((-(X + aa).^2)./c - (((Y + bb).^2)./20));

U1 = UL1 + UL2;
u = reshape(U1, 4096, 1);


end