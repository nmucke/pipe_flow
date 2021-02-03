clearvars
close all

%% inputs
% number of grid points
M  = 40;
% length of domain
Lx = 2*pi;

% note: filtering width = 2*H0
H0 = 2*Lx/M;

% number of points used to construct quadrature rule
nwidth = 2; % degree = 2*width+1
filter_type = 'cos';
%
t = 4*pi;

%% grid
dx = Lx/M;
x  = linspace(0,Lx,M+1)';
x  = x(1:end-1);

%% top-hat filtering of a sine signal

% exact signal:
u_ex = sin(x-t);

% filtered solution:
ubar  = -(cos(x+H0-t) - cos(x-H0-t))/(2*H0);

% nonuniform filter width, such that we have 'refinement' at the boundaries
% H     = (-0.5*cos(x)+1)*H0;
% derivative of non-uniform filter
% Hx    = (0.5*sin(x))*H0;
[H, Hx] = filter_width(x,H0,filter_type);

% non-uniform filtered solution:
ubar_nonuni = -(cos(x+H-t) - cos(x-H-t))./(2*H);

figure(1)
plot(x,u_ex,'s-');
hold on
plot(x,ubar,'o-')
grid
plot(x,ubar_nonuni,'x-')
legend('original','filtered','non-uniform filtered');

%% derivative of filtered quantities

% derivative of exact signal
dudx = cos(x-t);

% uniform filtering:
dubardx = (sin(x+H-t) - sin(x-H-t))./(2*H);
% note: bardu_dx = dubardx (differentiation and filtering commutes)

% non-uniform filtering
dubar_nonuni_dx = Hx./(2*H.^2) .*(cos(x+H-t) - cos(x-H-t)) + ...
    1./(2*H) .*((1+Hx).*sin(x+H-t) - (1-Hx).*sin(x-H-t));

% bar_dudx:
bar_nonuni_dudx = (sin(x+H-t) - sin(x-H-t))./(2*H);

% commutator error between dubar2dx and bar2dudx
% comm_error = dubar2dx - bar2dudx;
comm_error = -Hx.*ubar_nonuni./H + ...
    Hx./(2*H).*(sin(x+H-t) + sin(x-H-t));

% second derivative
% d2ubar_nonuni_dx2 = Hx./(2*H.^2) .*(cos(x+H-t) - cos(x-H-t)) + ...
%     1./(2*H) .*((1+Hx).*sin(x+H-t) - (1-Hx).*sin(x-H-t));

% commutator error based on Taylor expansion with u_xx
d2udx2 = -sin(x-t);
% comm_error_Taylor = (1/3).*Hx.*H.*d2udx2;


% commutator error based on Taylor expansion with bar{u}_xx
syms xx;
if (strcmp(filter_type,'cos'))
    g =@(xx) (-0.5*cos(xx)+1)*H0;
elseif (strcmp(filter_type,'uni'))
    g =@(xx) H0;
end

% filtered solution:
v   = @(xx) -(cos(xx+g(xx)-t) - cos(xx-g(xx)-t))./(2*g(xx));

% derivatives of filter and of filtered solution
gx  = diff(g,xx); % should match Hx
vx  = diff(v,xx); % should match dubar_nonuni_dx
vxx = diff(vx,xx); 

comm_error_Taylor = (1/3)*gx*g*vxx;

error_comm  = max(abs(double(subs(comm_error_Taylor,x))-comm_error))

figure(2)
plot(x,dudx,'s-')
hold on
plot(x,dubar_nonuni_dx,'o-');
plot(x,bar_nonuni_dudx,'x-');
legend('du/dx, unfiltered','d/dx u_{bar} (non-uniform filtered)','bar du/dx (non-uniform filtered)');

figure(3)
plot(x,comm_error,'s-');
hold on
plot(x,dubar_nonuni_dx-bar_nonuni_dudx,'x-')
% plot(x,comm_error_Taylor,'o-');
fplot(comm_error_Taylor,[0 Lx],'-');
legend('commutator error','commutator error', 'commutator error Taylor expansion');

figure(4)
plot(x,d2udx2);
hold on
fplot(vxx,[0 Lx]);
legend('u_{xx}, (unfiltered)','\bar{u}_{xx} (non-uniform filtered)')

%% create a quadrature rule with non-uniform filter width

addpath('/Users/sanderse/Dropbox/work/Programming/UQ/quadrature_rules');

% weight matrix
W     = zeros(M,M);

% extend spatial coordinates for periodicity
x_ext = [x-2*pi; x; x+2*pi];

for i=1:M

    % stencil
    stencil = i-nwidth:i+nwidth;
    % index in matrix
    indx    = mod(stencil,M);
    indx(indx==0) = M;
    
    % index in extended coordinates
    i_ext = i+M;
    indx2 = i_ext-nwidth:i_ext+nwidth;
    % point set used:
    x0    = x_ext(indx2); 
    
    % filter; note width = 2*D
    [H, Hx] = filter_width(x(i),H0,filter_type);
    
    % integration domain
    domain = [x(i)-H,x(i)+H];
    
    % row in weight matrix:
    W(i,indx) = QRgetWeights(length(x0),x0,domain);
    
end

figure(1)
ubar_quad = W*u_ex;

plot(x,ubar_quad,'d-');
legend('original','filtered','non-uniform filtered','non-uniform quadrature');