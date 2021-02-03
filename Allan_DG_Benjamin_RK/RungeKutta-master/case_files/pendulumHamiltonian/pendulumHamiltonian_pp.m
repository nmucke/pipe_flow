
plotopt = {'linewidth',1,'markersize',8};


theta = u_plot(1,:);
p     = u_plot(2,:);


%% solution
figure(1)
plot(t_plot,theta,'s-',plotopt{:})
legend(options.RK_method.name);
xlabel('t [s]');
ylabel('\theta [rad]');

%% phase portrait
figure(2)
plot(theta,p,'s-',plotopt{:})
xlabel('theta')
ylabel('p');


%% energy plots
figure(3)
dthetadt = p/(l^2);
V0 = l*g*(1-cos(theta(1)));
V  = l*g*(1-cos(theta)); % g*cos(theta); minus sign because highest potential energy for theta=pi
K0 = 0.5*(l*dthetadt(1))^2;
K  = 0.5*(l*dthetadt).^2; % 0.5*v^2, v=l*dtheta/dt

% in this case, Hamiltonian equals total energy

plot(t_plot,V-V0+K-K0)
legend('error in energy');



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
