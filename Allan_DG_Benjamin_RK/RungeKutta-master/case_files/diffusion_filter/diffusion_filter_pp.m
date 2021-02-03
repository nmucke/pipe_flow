%% diffusion PDE
% du/dt = nu*d2u/dx^2
% discretized as
% du/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

%% solution
figure(1)
surf(t_plot,x_in,u_plot)
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')


%%
figure(2)
plot(x_in,u_start,'s-');
hold on
plot(x_in,u_plot(:,end),'x-')
u_ex = diffusion_filter_ex(t_end,options);
plot(x_in,u_ex,'s-');
% plot(x_in,
legend('initial',options.RK_method.name,'exact');

err(j) = max(abs(u_ex - u_plot(:,end)));
%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on