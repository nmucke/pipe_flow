%% showme.me
%
% Show data for shallow water problem.

figure
subplot(1,4,1)
plot(x0,h0,'b')
xlabel('Length along channel, x')
ylabel('Initial condition, free surface elevation')
axis([min(x0(:)) max(x0(:)) -0.1 max(h0(:))*1.1])
grid on

subplot(1,4,2)
plot(x,h)
xlabel('Length along channel, x')
ylabel('Height, h')
axis([min(x(:)) max(x(:)) -0.1 max(h(:))*1.1])
grid on

subplot(1,4,3)
plot(x,u)
xlabel('Length along channel, x')
ylabel('Velocity, u')
axis([min(x(:)) max(x(:)) -0.1 max(u(:))*1.1])
grid on

subplot(1,4,4)
plot(x,hu)
xlabel('Length along channel, x')
ylabel('Momentum, u*h')
axis([min(x(:)) max(x(:)) -0.1 max(hu(:))*1.1])
grid on

set(gcf,'position',[120 532 1726  513])
