
figure(2)
plot(u_plot(1,:),u_plot(2,:),'s-')
hold on

figure(3)
% get total energy and angular momentum
H0 = -1/2;
L0 = sqrt(1-e^2);
H  = 0.5*(u_plot(3,:).^2 + u_plot(4,:).^2) - 1./sqrt(u_plot(1,:).^2 + u_plot(2,:).^2);
L  = u_plot(1,:).*u_plot(4,:) - u_plot(2,:).*u_plot(3,:);
plot(t_plot,abs(H-H0),'s-');
hold on
plot(t_plot,abs(L-L0),'o-');
legend('error in energy','error in angular momentum');


