%% linear ODE
% du/dt = lambda*u
% exact solution u(t) = u(t0)*exp(lambda*t);

%% solution
figure(1)
plot(t_plot,u_plot,'s-')
hold on
plot(t_plot,simple_ex(t_plot,options),'o-')

legend(options.RK_method.name,'exact');

%% error plots

figure(2)

for k=1:RK_number
loglog(1./N_list,err(:,k),'s-');
hold on
get_slope(1./N_list,err(:,k),Nsim)
end
grid on
legend(RK_list)
