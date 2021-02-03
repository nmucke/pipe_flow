function u = diffusion_filter_ex(t,options)
% du/dt = nu*d2u/dx^2

% alpha = options.model.constants(1); % diffusion coefficient
% L     = options.model.constants(2); % domain length
x_in  = options.model.x_in;
u     = t + sin(2*pi*x_in) + sin(8*pi*x_in);


end