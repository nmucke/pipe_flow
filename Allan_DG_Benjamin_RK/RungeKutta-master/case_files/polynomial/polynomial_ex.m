function u = polynomial_ex(t,options)
% polynomial in time
% du/dt = t^p
% exact solution u(t) = u(t0) + t^(p+1)/(p+1);

p       = options.model.constants(1);
u0      = options.model.u0;
t_start = options.time.t_start;

u       = u0 + (t.^(p+1) - t_start.^(p+1))/(p+1); 

end