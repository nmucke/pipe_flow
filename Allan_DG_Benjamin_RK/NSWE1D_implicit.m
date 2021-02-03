function [h,hu] = NSWE1D_implicit(h,hu,FinalTime,g,dt)
%
% function [h,hu] = NSWE1D(h,hu,FinalTime,g)
% Purpose : Integrate 1D NSWE until FinalTime starting with
%           initial the condition h and hu
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
Globals1D;

time = 0;

un = [h,hu];
un = un(:);

[t,y] = ode45(@(t,u) NSWERHS_implicit(u,t),[0,2.5],un);


% 
% % outer time step loop
% while time < FinalTime
% 
%     un = RKstep_dirk(un,time,dt);
%     
%     
%     % Increment time
%     time = time+dt;
%     
%     h = un(1:end/2);
%     hu = un(end/2+1:end);
%     
%     
%     
%     % Visualization
%     subplot(2,1,1)
%     plot(x(:),h(:))
%     title(sprintf('Solution at t = %.2f',time))
%     xlabel('x') , ylabel('h')
%     %ylim([1 4])
% 
%     subplot(2,1,2)
%     plot(x(:),hu(:))
%     xlabel('x') , ylabel('hu')
%     %ylim([-5 10])
%     drawnow
%     %pause
% end

h = y(:,1:end/2);
hu = y(:,end/2+1:end);


% Visualization
subplot(2,1,1)
plot(x(:),h(end,:))
title(sprintf('Solution at t = %.2f',time))
xlabel('x') , ylabel('h')
%ylim([1 4])

subplot(2,1,2)
plot(x(:),hu(end,:))
xlabel('x') , ylabel('hu')
%ylim([-5 10])
drawnow
%pause

return

