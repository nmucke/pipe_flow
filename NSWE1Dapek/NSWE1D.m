function [h,hu] = NSWE1D(h,hu,FinalTime,g)
%
% function [h,hu] = NSWE1D(h,hu,FinalTime,g)
% Purpose : Integrate 1D NSWE until FinalTime starting with
%           initial the condition h and hu
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
Globals1D;

time = 0;

% Runge?Kutta residual storage
resh  = zeros(Np,K);
reshu = zeros(Np,K);

% Compute time step size
u = hu./h;

lambda = max(max(abs([u+sqrt(g*h/2) ; u-sqrt(g*h/2)])));
CFL = 0.2;
xmin = min(abs(x(1,:)-x(2,:)));
dt = CFL*xmin/lambda;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps

    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsh,rhshu] = NSWERHS1D(h,hu,timelocal,g);
        resh = rk4a(INTRK)*resh + dt*rhsh;
        reshu = rk4a(INTRK)*reshu + dt*rhshu;
        h = h+rk4b(INTRK)*resh;
        hu = hu+rk4b(INTRK)*reshu;
    end
    
    % Increment time
    time = time+dt;
    
    % Visualization
    subplot(2,1,1)
    plot(x(:),h(:))
    title(sprintf('Solution at t = %.2f',time))
    xlabel('x') , ylabel('h')
    ylim([1 4])

    subplot(2,1,2)
    plot(x(:),hu(:))
    xlabel('x') , ylabel('hu')
    ylim([-5 10])
    drawnow
    %pause
end
return

