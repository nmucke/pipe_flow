%% linear ODE
% du/dt = lambda*u
% exact solution u(t) = u(t0)*exp(lambda*t);

plotopt = {'linewidth',1,'markersize',8};

theta    = u_plot(1,:);
dthetadt = u_plot(2,:);

%% solution
figure(1)
plot(t_plot,theta,'s-',plotopt{:})
legend(options.RK_method.name);
xlabeltex('t [s]',14);
ylabeltex('\theta [rad]',14);
grid

%% phase portrait
figure(2)
plot(theta,dthetadt,'s-',plotopt{:})
xlabeltex('\theta',14);
ylabeltex('\dot{\theta}',14);
grid

%% energy plots
figure(3)
% V0 = l*g*(1-cos(u_start(1)));
V  = l*g*(1-cos(theta)); % g*cos(theta); minus sign because highest potential energy for theta=pi
% K0 = 0.5*(l*u_start(2))^2;
K  = 0.5*(l*dthetadt).^2; % 0.5*v^2, v=l*dtheta/dt

plot(t_plot,V+K-(V(1)+K(1)))
legend('error in energy');

%%
plot(t_plot,V+K,plotopt{:})
xlabeltex('t',14)
ylabeltex('E',14)
grid

%% error plots

% figure(4)
% 
% for k=1:RK_number
% loglog(1./N_list,err(:,k),'s-');
% hold on
% get_slope(1./N_list,err(:,k),Nsim)
% end
% grid on
% legend(RK_list)
