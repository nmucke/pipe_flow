%% PDE:
% du/dt + u du/dx = nu*d2u/dx^2
% or equivalently
% du/dt + d (0.5*u^2) /dx = nu*d2u/dx^2
% discretized as
% du/dt + C*(u^2) = D * u


%% solution
figure(1)
surf(t_plot,x_in,u_plot)
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')


%%
figure(2)
plot(x_in,u_plot(:,end),'s-')
hold on
u_ex = Burgers_ex(t_end,options);
plot(x_in,u_ex);
% plot(x_in,
legend(options.RK_method.name,'exact');


%% error plots

if (Nsim>1)
    figure(3)
    loglog(1./N_list,err,'s-');
    get_slope(1./N_list,err,Nsim)
    grid on
end