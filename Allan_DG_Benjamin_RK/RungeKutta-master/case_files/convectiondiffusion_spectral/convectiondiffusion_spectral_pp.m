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
plot(x_in,u_plot(:,end),'x-')
hold on
u_ex = convectiondiffusion_spectral_ex(t_end,options);
plot(x_in,u_ex,'s-');
% plot(x_in,
legend(options.RK_method.name,'exact');


%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on