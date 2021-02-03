function u = Burgers_filter_ex(t,options,x)
% BURGERS_EX exact solution to Burgers equation
% du/dt + d/dx(0.5* u^2) = nu*d2u/dx^2
% see Rang, JCAM 2015, Improved traditional Rosenbrock–Wanner methods for stiff
% ODEs and DAEs; also Jens Lang; Jan Verwer; Whitham ch. 4

% argument x is optional; if not provided, the solution at all interior
% points x is given

nu = options.model.constants(1); % diffusion coefficient

if (nargin<3)
    x  = options.model.x;
end

% Rang / Lang / Verwer solution:
r1 = exp(-(x-0.5)/(20*nu) - 99*t/(400*nu));
r2 = exp(-(x-0.5)/(4*nu) - 3*t/(16*nu));
r3 = exp(-(x-0.375)/(2*nu));

% alternative, own solution, based on book Whitham, ch. 4:
% c1 = 1/10;
% c2 = 1/2;
% c3 = 1;
% 
% b1 = -1/2*c1/(2*nu);
% b2 = -1/2*c2/(2*nu);
% b3 = -3/8*c3/(2*nu);
% 
% r1 = exp(-c1*x/(2*nu) + (c1^2)*t/(4*nu) - b1);
% r2 = exp(-c2*x/(2*nu) + (c2^2)*t/(4*nu) - b2);
% r3 = exp(-c3*x/(2*nu) + (c3^2)*t/(4*nu) - b3);

u  = (0.1*r1 + 0.5*r2 + r3) ./ (r1 + r2 + r3);

%% check the exact solution:
% syms x t v(x,t) u(x,t)
% nu = 0.01;
% c1 = 1/10;
% c2 = 1/2;
% c3 = 1;
% 
% b1 = -1/2*c1/(2*nu);
% b2 = -1/2*c2/(2*nu);
% b3 = -3/8*c3/(2*nu);
% 
% r1 = exp(-c1*x/(2*nu) + (c1^2)*t/(4*nu) - b1);
% r2 = exp(-c2*x/(2*nu) + (c2^2)*t/(4*nu) - b2);
% r3 = exp(-c3*x/(2*nu) + (c3^2)*t/(4*nu) - b3);

% In Rang / Lang / Verwer, the following is used:
% r1 = exp(-(x-1/2)/(20*nu) - 99*t/(400*nu));
% r2 = exp(-(x-1/2)/(4*nu) - 3*t/(16*nu));
% r3 = exp(-(x-3/8)/(2*nu));
% 
% u(x,t)  = ( (1/10)*r1 + (1/2)*r2 + r3) ./ (r1 + r2 + r3);
% 
% v(x,t)  = simplify(diff(u,t) + u*diff(u,x) - nu*diff(diff(u,x),x));

% % t_test = 0:0.01:1;
% x_test = 0:0.01:1;
% t_test = 0:0.1:1;
% for i=1:length(t_test)
%     plot(x_test,u(x_test,t_test(i)))
% %     plot(x_test,v(x_test,t_test(i)))
%     hold on
% end


end