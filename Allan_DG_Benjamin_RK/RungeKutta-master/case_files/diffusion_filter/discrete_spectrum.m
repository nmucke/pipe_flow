function [kx,G_abs,G_real,G_imag] = discrete_spectrum(coeffs,offset,Mk,plot_output)
% get amplification factor (spectral analysis) for quadrature rules
% examples:
% discrete_spectrum([1/6; 2/3; 1/6],[-1;0;1],100);
% or (see thesis Edoh, figure 2.3):
% [kx,Gabs,Greal,Gimag]=discrete_spectrum([1/16; 3/4; 3/8; -1/4; 1/16],[-1;0;1;2;3],100,1);

% coeffs  = [1/4; 1/2; 1/4]; % trapz
% coeffs  = [1/6; 2/3; 1/6]; % simpson
% offset  = [-1; 0; 1]; % stencil of points


if (nargin<2)
    error('too few arguments');
end

% Mk: number of points
if (nargin<3)
    Mk = 100;
end

if (nargin<4)
    plot_output = false;
end

% discretization of k*deltax on [0,1] with Mk points
kx = linspace(0,2*pi,Mk)';

j  = sqrt(-1);

G_abs = zeros(Mk,1);
G_real = zeros(Mk,1);
G_imag = zeros(Mk,1);

for k=1:Mk
    
    G        = sum(coeffs.*exp(kx(k)*offset*j));
    G_abs(k) = abs(G);
    G_real(k) = real(G);
    G_imag(k) = imag(G);
    
end

if (plot_output)
    figure
    plot(kx/pi,G_abs,'s-');
    hold on
    plot(kx/pi,G_real,'o-');
    plot(kx/pi,G_imag,'d-');
        plot(kx/pi,G_exact,'s-');

    legend('absolute value','real part','imaginary part','exact')
end