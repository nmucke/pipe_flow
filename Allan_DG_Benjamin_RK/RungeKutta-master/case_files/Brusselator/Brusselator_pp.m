%%  PDE
% du/dt = A + u^2*v - (B+1)*u + alpha*d2u/dx2
% dv/dt = B*u - u^2*v + alpha*d2v/dx2

% discretized as
% du/dt = A + u^2*v - (B+1)*u + alpha*D*u + uBC
% dv/dt = B*u - u^2*v + alpha*D*v + vBC

%% solution
figure(1)
surf(t_plot,x_in,u_plot(options.model.ind_u,:))
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')


%%
% figure(2)
% plot(x_in,u_plot(:,end))
% hold on
% u_ex = diffusion_ex(t_end,options);
% plot(x_in,u_ex);
% % plot(x_in,
% legend(options.RK_method.name,'exact');


%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on