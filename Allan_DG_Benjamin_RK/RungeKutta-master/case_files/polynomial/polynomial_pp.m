%% polynomial in time
% du/dt = t^p
% exact solution u(t) = u(t0)*t^(p+1)/(p+1);

%% solution
figure(1)
plot(t_plot,u_plot,'s-')
hold on
plot(t_plot,polynomial_ex(t_plot,options),'o-')

legend(options.RK_method.name,'exact');



%% error plots

figure(2)
loglog(1./N_list,err,'s-');
get_slope(1./N_list,err,Nsim)
grid on