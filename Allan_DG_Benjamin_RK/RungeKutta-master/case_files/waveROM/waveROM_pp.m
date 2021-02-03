%% wave PDE
% d2u/dt = c*d2u/dx^2
fontsize = 14;
fontname = 'Helvetica';
set(0,'defaultlinelinewidth',2)


if (strcmp(ROM_basis,'symplectic'))
    U_FOM    = Aplus'*u_plot; % [q(...); p(...)]
elseif (strcmp(ROM_basis,'POD'))
    U_FOM    = V*u_plot;
end

q = U_FOM(1:Morig,:);
p = U_FOM(Morig+1:end,:);

q_exact = snapshots.q(:,1:N+1);

% theta1  = U_FOM(1,:);
% theta2  = U_FOM(2,:);
% p1      = U_FOM(3,:);
% p2      = U_FOM(4,:);


%% solution
figure(1)
set(gcf,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
surf(t_plot,x_in,q)
shading interp
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')
zlabel('u')
title(['ROM ' ROM_basis])


%%
figure(2)
set(gcf,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
plot(x_in,q_exact(:,N+1))
hold on
plot(x_in,q(:,end))
% u_ex = diffusion_ex(t_end,options);
% plot(x_in,u_ex);
% plot(x_in,
% legend(options.RK_method.name,'exact');
axis([0 1 0 1])
grid
xlabel('x')
ylabel('u')


%% solution accuracy
q_error = sqrt((1/Morig)*sum((q-q_exact).^2));
figure
set(gcf,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
semilogy(t_plot,q_error)
xlabel('t')
ylabel('solution error')
grid


%% Hamiltonian accuracy

% 'exact' (there is still a discretization error) H from high fidelity
% dataset
H_exact = snapshots.H;

% discrete Hamiltonian
% circshift shifts along columns
qnext = circshift(q,-1);
qprev = circshift(q,1);
% sum operates over the first (=spatial) dimension
H     = sum(dx*(0.5*p.^2 + (c^2)/(4*dx^2) * ( (qnext - q).^2 + (q - qprev).^2 ) ) );

figure
set(gcf,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
plot(t_plot,H_exact(1:N+1));
hold on
plot(t_plot,H); %-H_exact(1));

grid
xlabel('t')
ylabel('H')

%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on