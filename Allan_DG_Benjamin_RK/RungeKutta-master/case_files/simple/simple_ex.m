function u = simple_ex(t,options)
% linear ODE
% du/dt = lambda*u
% exact solution u(t) = u(t0)*exp(lambda*t);

lambda  = options.model.constants(1);
u0      = options.model.u0;
t_start = options.time.t_start;

u       = u0*exp(lambda*(t-t_start)); 
