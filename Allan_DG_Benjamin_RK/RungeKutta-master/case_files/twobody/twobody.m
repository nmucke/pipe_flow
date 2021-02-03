% twobody:
M       = 4; % dimension of problem 
options.model.function = 'twobody';
e = 0.6; % eccentricity
u_start = [1-e;0;0;sqrt((1+e)/(1-e))]; % initial condition: posx,posy,velx,vely

N_list = 200;
t_start = 0;
t_end   = 20;

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 2;