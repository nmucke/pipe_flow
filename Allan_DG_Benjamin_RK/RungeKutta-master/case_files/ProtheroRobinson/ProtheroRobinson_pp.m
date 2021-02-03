%% Prothero-Robinson problem
% Prothero & Robinson,﻿On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations
% the original problem is
% du/dt = g'(t) + lambda*(u - g(t))
% exact solution  u = g(t)

% here we take a slightly different problem:
% du/dt = -lambda*(u - cos(t))
% exact solution u(t) = (lambda^2/(lambda^2+1))*(-exp(-lambda*t) + cos(t) + sin(t)/lambda)


%% solution
figure(1)
plot(t_plot,u_plot,'s-')
hold on
plot(t_plot,ProtheroRobinson_ex(t_plot,options),'o-')

legend(options.RK_method.name,'exact');



%% error plots
if (Nsim>1)
    figure(2)

    for k=1:RK_number
    loglog(1./N_list,err(:,k),'s-');
    hold on
    get_slope(1./N_list,err(:,k),Nsim)
    end
    grid on
    legend(RK_list)
%     figure(2)
%     loglog(1./N_list,err,'s-');
%     get_slope(1./N_list,err,Nsim)
%     grid on
end