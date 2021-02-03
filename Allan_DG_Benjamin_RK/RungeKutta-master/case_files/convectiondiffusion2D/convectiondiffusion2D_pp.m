%% PDE:
% du/dt = -c*du/dx + alpha*d2u/dx^2
% discretized as
% du/dt = A*u = -C*u + D*u
% with periodic boundary conditions

if (average==1)
    u_rec = Rec*u_plot;
    u_avg = u_plot;
else
    u_rec = u_plot;
    u_avg = Ave*u_plot;
end

%% (reconstructed) solution in space at different time steps

if (options.pp.movie == 1)
    writerObj = VideoWriter(options.pp.moviename,'Motion JPEG AVI');
    writerObj.FrameRate = options.pp.movierate;
    open(writerObj);
end

figure
for i=1:floor(N/50):N
    
    title('reconstructed solution');
    surf(x_in,y_in,reshape(u_rec(:,i),Mx,My)');
    zlim([-1 3])
    xlabel('x')
    ylabel('y')
    zlabel('u')
    
    if (options.pp.movie == 1)
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    
    %     pause(0.1)
    %     pause;
    
end

if (options.pp.movie == 1)
    close(writerObj);
end
% figure
%
% set(gcf,'color','w');
% surf(x_in,t_plot,u_plot')
% % surf(t_plot,x_in,u_plot)
% % hold on
% % plot(t_plot,polynomial_ex(t_plot,options),'o-')
% xlabel('x')
% ylabel('t')
% zlabel('u')
% colorbar
%
% set(gca,'LineWidth',1)
% set(gca,'FontSize',13)

%%
ix = Mx;
u_rec_3D = reshape(u_rec,Mx,My,N+1);
% get solution at ix location
u_rec_2D = squeeze(u_rec_3D(ix,:,:));

iy = ceil(My/2);
if (iy == My/2)
    u_rec_1D = (u_rec_2D(iy,:) + u_rec_2D(iy+1,:))/2;
else
    u_rec_1D = u_rec_2D(iy,:);
end

figure
% open(
plot(t_plot,u_rec_1D)
xlabel('t')
ylabel('reconstructed solution at x=1, y=1/2')
grid

%% averaged solution at outlet
if (average==0)
    My_new = 1;
end
% My_new = size(Av1D,1);
u_avg_3D    = reshape(u_avg,Mx,My_new,N+1); % reshape to 3D array
u_avg_2D    = squeeze(sum(u_avg_3D,2)/My_new); % take sum in y-direction to sum averages at each x-location (to get approximation to int u dy from 0..1)

figure
plot(t_plot,u_avg_2D(Mx,:))
xlabel('t')
ylabel('average solution at x=1')
grid

%%
% figure
% plot(x_in,u_plot(:,end),'x-')
% hold on
% u_ex = convectiondiffusion_ex(t_end,options);
% plot(x_in,u_ex,'s-');
% % plot(x_in,
% legend(options.RK_method.name,'exact');


%% error plots

% figure(2)
% loglog(1./N_list,err,'s-');
% get_slope(1./N_list,err,Nsim)
% grid on