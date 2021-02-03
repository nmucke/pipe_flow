%% diffusion PDE
% du/dt = nu*d2u/dx^2
% discretized as
% du/dt = D * u
% with zero Dirichlet boundary conditions at x=0 and x=L

% u_plot is the ROM solution
% map to FOM space
u_FOM    = V*u_plot;

u_exact = snapshots.u_plot(:,1:N+1);


%% solution
figure
surf(t_plot,x_in,u_FOM)
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('t')
ylabel('x')


%%
figure
plot(x_in,u_FOM(:,end),'x-')
hold on
u_ex = snapshots.u_ex; % diffusionROM_ex(t_end,options);
plot(x_in,u_exact(:,end),'s-');
plot(x_in,u_ex,'s-');
% plot(x_in,
legend('ROM','FOM','exact');


%%
figure
x_plot = options.model.Morig/2;
plot(t_plot,u_FOM(x_plot,:),'x-')
hold on
u_ex = snapshots.u_ex; % diffusionROM_ex(t_end,options);
plot(t_plot,u_exact(x_plot,:),'s-');
xlabel('t')
ylabel('u');
legend('ROM','FOM');


%% error analysis
% error in space and time
e = u_FOM - u_exact;
max(max(abs(e)))