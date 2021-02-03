clearvars
close all

Mk = 100;
plot_output = 0;

fg_ratio = 2; % Delta / dx;

% top-hat filter with Delta
filter = 'box';
% get G(k*Delta)
[kD_tophat, G_tophat] = continuous_spectrum(filter,Mk,plot_output);

% Gaussian filter with Delta
filter = 'gauss';
% get G(k*Delta)
[kD_gauss, G_gauss] = continuous_spectrum(filter,Mk,plot_output);

% Simpson
coeff   = [1/6; 2/3; 1/6];
stencil = [-1;0;1]; % corresponding to -Delta x, 0, Delta x
% we get G(k*Deltax) but want to plot G(k*Delta) = G(2*k*Deltax)
[k_simps, G_simps] = discrete_spectrum(coeff,stencil,Mk,plot_output);

% trapz
coeff   = [1/4; 1/2; 1/4];
stencil = [-1;0;1]; % corresponding to -Delta x, 0, Delta x
% we get G(k*Deltax) but want to plot G(k*Delta) = G(2*k*Deltax)
[k_trapz, G_trapz] = discrete_spectrum(coeff,stencil,Mk,plot_output);



figure
plot(kD_tophat/pi, G_tophat,'s-','LineWidth',2)
hold on
plot(kD_gauss/pi, G_gauss,'s-','LineWidth',2)
plot(fg_ratio*k_simps/pi, G_simps,'o-','LineWidth',2);
plot(fg_ratio*k_trapz/pi, G_trapz,'o-','LineWidth',2);
xlim([0 2])
legend('Top-hat (continuous)','Gauss (continuous)','Simpson (discrete)','Trapezoidal (discrete)');
xlabel('k*Delta/pi');
ylabel('G');
set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)
grid