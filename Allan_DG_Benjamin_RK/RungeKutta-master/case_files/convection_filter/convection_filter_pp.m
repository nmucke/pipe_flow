%% PDE:
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 




%% solution
figure

set(gcf,'color','w');
surf(x_in,t_plot,u_plot')
% surf(t_plot,x_in,u_plot)
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('x')
ylabel('t')
zlabel('u')
colorbar

set(gca,'LineWidth',1)
set(gca,'FontSize',13)


%%
figure
plot(x_in,u_plot(:,end),'s-')
hold on
u_ex = convection_filter_ex(t_end,options);
% u_ex = sin(x_in - t_end);
plot(x_in,u_ex,'x-');
% plot(x_in,
legend('N=20','exact filtered');
grid
xlabel('x')
ylabel('u')

max(abs(u_plot(:,end)-u_ex))

%% error plots
% err = max(abs(u_ex - unew))
% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on

%%
% N
% 20 40 80
% error saturates because quadrature rule is over smaller domain
%    0.092431491942225; 0.023233210297442 0.014470840960025

% 3-point stencil (D=2)
% N
% 20 40 80
% error
%     0.107841807052903  0.086666584834545 0.082711627776459

% varying point stencil (skip=2)
% 
% N = [20 40 80];
% error = [0.107841807052903  0.023232226673180 0.005814617505442];
% loglog(N, error,'s-')
% ylabel('max. error at t=1')
% xlabel('N');
% xlim([10 100])
% grid

%% comparison with classic approach
% fixed filter width L/4
% t_end = pi/32
% N = [20 40 80];
% error_classic = [0.029209648622271 0.031065342235109  0.032007749134922];
% skip=2:
% error_new = [0.099686719763016  0.010599248597383  9.076501949475957e-05];
% skip=1:
% error_new = [0.011670163385267 3.614377401734217e-04 diverges]

% figure
% loglog(N, error_new,'s-')
% hold on
% loglog(N, error_classic,'o-')
% ylabel('max. error at t=pi/32')
% xlabel('N');
% legend('new approach','classic approach');
% grid

%% comparison with classic approach
% varying filter width: dx
% 3-point quadrature:
error_new = [0.102246223854682 0.025778163327323  0.006456755409862];
error_classic = [0.099723763204308 0.025659181777464  0.006449617342576];

figure
loglog(N, error_new,'s-')
hold on
loglog(N, error_classic,'o-')
ylabel('max. error at t=2*pi')
xlabel('N');
legend('new approach','classic approach');
xlim([10 100])
grid



