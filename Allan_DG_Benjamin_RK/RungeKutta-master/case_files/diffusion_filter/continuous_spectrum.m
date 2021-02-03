function [kD,G_abs] = continuous_spectrum(filter,Mk,plot_output)
% get amplification factor (spectral analysis) for continuous filters
% examples:
% filter.type = 'gauss'; filter.gamma = 6;
% continuous_spectrum(filter,100,1);
% or 
% filter='spectral';
% continuous_spectrum(filter,100,1);


if (nargin<1)
    error('too few arguments');
end
if (isstruct(filter))
    filter_type = filter.type;
else
    filter_type = filter;
end

% Mk: number of points
if (nargin<2)
    Mk = 100;
end

if (nargin<3)
    plot_output = false;
end

% discretization of k*Delta on [0,pi] with Mk points
kD = linspace(0,2*pi,Mk)';

G     = zeros(Mk,1);
G_abs = zeros(Mk,1);
% G_real = zeros(Mk,1);
% G_imag = zeros(Mk,1);


for k=1:Mk
    
    switch filter_type
        
        case 'gauss'
            % kernel:
            % G =@(r) sqrt( (gamma/(pi* (Delta^2)))) * exp(- gamma*(r.^2) / Delta^2);
            % %
            % int( exp(k*ksi*i) * G(x-ksi), ksi, -inf, inf) =
            %  = sqrt(a/pi)*exp(-a*x^2)*int(exp(-a*ksi^2 + (2*a*x+b*i)*ksi),ksi,-inf,inf)
            %  = sqrt(a/pi)*exp(-a*x^2)*(pi^(1/2)*exp((2*a*x + b*1i)^2/(4*a)))/a^(1/2)
            %  = exp(-a*x^2)*exp((2*a*x + b*1i)^2/(4*a))
            % where a = gamma/Delta^2, b = k
            % take x = 0, then
            % = exp(-b^2/(4*a)) = exp( - (k Delta)^2 / (4*gamma))
            if (isfield(filter,'gamma'))
                gamma = filter.gamma;
            else
                gamma = 6;
            end
            %             Delta = filter.Delta;
            G = exp(- kD(k)^2/(4*gamma));
            
        case 'spectral'
            %  G =@(r) (sin(pi*(r+eps)/Delta)./(pi*(r+eps)));
            G = kD(k)<=pi;
            
        case 'box'
            % kernel:
            % G =@(r) (1/Delta).*(abs(r)<=Delta/2);
            % (1/Delta)*int( exp(k*ksi*i) * G(x-ksi), ksi, x-D/2, x+D/2) =
            % (1/Delta)*(1/(i*k))*(exp(k*Delta*i/2) - exp(-k*Delta*i/2)) = 
            % sin(k*Delta/2)/(k*Delta/2)
            G = sin(-kD(k)/2)/(kD(k)/2);
            
    end
    
    G_abs(k) = abs(G);
    %     G_real(k) = real(G);
    %     G_imag(k) = imag(G);
    
end

if (plot_output)
    figure
    plot(kD,G_abs,'s-');

end