%% wave PDE
% d2u/dt = c*d2u/dx^2

plotopt = {'linewidth',2,'markersize',12};

q = u_plot(1:M/2,:);
p = u_plot(M/2+1:end,:);

%% solution
figure(1)
surf(t_plot,x_in,q,plotopt{:})
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')
zlabel('u')

shading interp

%%
ind_t = 1+(2.5/dt);

figure(2)
% plot(x_in,q(:,end))
plot(x_in,q(:,ind_t))

% hold on
% u_ex = diffusion_ex(t_end,options);
% plot(x_in,u_ex);
% plot(x_in,
legend(options.RK_method.name); %,'exact');
axis([0 1 0 1])


%% energy / Hamiltonian
% discrete Hamiltonian
% circshift shifts along columns
qnext = circshift(q,-1);
qprev = circshift(q,1);
% sum operates over the first (=spatial) dimension
H     = sum(dx*(0.5*p.^2 + (c^2)/(4*dx^2) * ( (qnext - q).^2 + (q - qprev).^2 ) ) );

figure
plot(t_plot,H-H(1));
grid
xlabel('t')
ylabel('error in H')

%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on