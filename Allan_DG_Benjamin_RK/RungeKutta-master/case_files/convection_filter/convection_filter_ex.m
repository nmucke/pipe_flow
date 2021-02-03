function u = convection_filter_ex(t,options)
% du/dt = -c*du/dx 
% discretized as
% du/dt = A*u = -C*u 
% with periodic boundary conditions 
% and IC sin(x)

c     = options.model.constants(1); % convection coefficient
H0    = options.model.constants(2); % filter width
k     = options.model.constants(3); % convection coefficient

filter_type = options.model.filter_type; % filter type
x_in  = options.model.x_in;
% H0    = options.model.H0;

% exact solution:
% u_ex  = sin((x_in-c*t));

% filtered solution:
[H, ~] = filter_width(x_in,H0,filter_type);
u      = -(cos(k*(x_in+H)-c*t) - cos(k*(x_in-H)-c*t))./(2*k*H);

% u_start = -(cos(k*(x_in+H)-t_start) - cos(k*(x_in-H)-t_start))./(2*k*H);


end