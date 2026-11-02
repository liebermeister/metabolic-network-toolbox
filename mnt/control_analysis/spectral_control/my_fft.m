function [omega,y] = my_fft(x,T)

% let x contain the values of x(t) at
% t = -T .. 0 .. T
%
% calculate values of y and omega such that
% y(omega) is (approximately) the continuous Fourier transform of x(t)

if ~exist('T','var'), T=1; end

dt = 2*T / (length(x)-1 );
y  = fftshift(fft(ifftshift(x)));

y     = y * dt;
f     = -1/dt : 1/T : 1/dt;
omega = f/2*(2*pi);
