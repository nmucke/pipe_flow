%% PDE:
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions 

% u_plot is the ROM solution
% map to FOM space
u_FOM    = V*u_plot;
% u_exact is the original full order model solution
u_exact = snapshots.u_plot(:,1:N+1);
% u_ex is the exact solution to the PDE (no discretization error):
u_ex = snapshots.u_ex; 


%% full space-time solution
figure
surf(t_plot,x_in,u_FOM)
% hold on
% plot(t_plot,polynomial_ex(t_plot,options),'o-')
xlabel('x')
ylabel('t')


%% solution at end time
figure
plot(x_in,u_FOM(:,end),'x-')
hold on
plot(x_in,u_exact(:,end),'s-');
plot(x_in,u_ex,'s-');
legend('ROM','FOM','exact');


%% solution at a single point in space as a function of time
figure
x_plot = 1; %options.model.Morig/2;
plot(t_plot,u_FOM(x_plot,:),'x-')
hold on
u_ex = snapshots.u_ex; 
plot(t_plot,u_exact(x_plot,:),'s-');
xlabel('t')
ylabel('u');
legend('ROM','FOM');


%% error analysis
% error in space and time
e = u_FOM - u_exact;
max(max(abs(e)))
% error at final time step
e = u_FOM(:,end) - u_exact(:,end);
max(abs(e))