%% polynomial in time
% du/dt = t^p
% exact solution u(t) = u(t0)*t^(p+1)/(p+1);

%% solution
figure(1)
plot(t_plot,u_plot,'s-')
hold on
plot(t_plot,Trefethen1_ex(t_plot,options),'o-')

legend(options.RK_method.name,'exact');



%% error plots

figure(2)
for k=1:RK_number
loglog(1./N_list,err(:,k),'s-');
hold on
disp(RK_list{k});
get_slope(1./N_list,err(:,k),Nsim)
end
grid on
legend(RK_list)

% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on