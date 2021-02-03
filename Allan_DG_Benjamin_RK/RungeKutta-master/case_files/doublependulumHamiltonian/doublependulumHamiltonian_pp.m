%% linear ODE
% du/dt = lambda*u
% exact solution u(t) = u(t0)*exp(lambda*t);

theta1  = u_plot(1,:);
theta2  = u_plot(2,:);
p1  = u_plot(3,:);
p2  = u_plot(4,:);


plotopt = {'linewidth',2,'markersize',12};

%% solution
figure(1)

x1 = L1*sin(theta1);
y1 = -L1*cos(theta1);
x2 = x1 + L2*sin(theta2);
y2 = y1 - L2*cos(theta2);

plot(x1,y1,'-',plotopt{:})
hold on
plot(x2,y2,'-',plotopt{:})
xlim([-(L1+L2) (L1+L2)]);
ylim([-(L1+L2) (L1+L2)]);
axis square
grid
% plot(t_plot,u_plot(1,:),'s-')
% legend(options.RK_method.name);
% xlabel('t [s]');
% ylabel('theta [rad]');

%% phase portrait
figure(2)
plot(u_plot(1,:),u_plot(2,:),'s-')
xlabel('\theta_1')
ylabel('\theta_2')
axis square
grid

%% energy plots
figure(3)
V = -(m1+m2)*g*L1*cos(theta1) - m2*g*L2*cos(theta2);
K = ((L2^2)*m2*(p1.^2) + (L1^2)*(m1+m2)*(p2.^2) - 2*m2*L1*L2*p1.*p2.*cos(theta1-theta2))./...
    (2*(L1^2)*(L2^2)*m2*(m1+m2*sin(theta1-theta2).^2));
plot(t_plot,V+K-(V(1)+K(1)))
legend('Error in Hamiltonian');

