%% RK_generic
% time integration of systems of equations du/dt = f(u(t),t)
% with generic (explicit, implicit, dirk) Runge-Kutta methods

% steps to solve an ODE with RK_generic:
% 0. decide on a project name or use existing project from case_files
% directory
% 1. create a folder in the case_files directory that has as name the project name
% 2. create at least the following two files: project.m and project_F.m,
% where project stands for your project name. A simple way of starting is
% by copying an existing folder and changing names. In project.m you
% specify parameters such as start time, end time, time step, etc.; in
% project_F.m you specify the right-hand side f(u(t),t) of the ODE du/dt=f
% 3. optionally, create the files project_ex.m (if you have an exact or
% reference solution) and project_pp.m (if you want to do postprocessing after
% simulation is finished)
% 4. set case_name (couple of lines below here) to your project name
% 5. choose the RK method or list of RK methods below
% 6. run this file

clear all
close all
clc


global rhs_evals;


% if jacobian==2, automatic differentiation package is needed
autodiff_type = 2; % 1: AutoDiff_R2016b, 2: automatic differentiation (preferred)
path_adi_1 = '/Users/sanderse/Dropbox/work/Programming/libs/AutoDiff_R2016b/';
path_adi_2 = 'automatic_differentiation';

warning('off','MATLAB:rmpath:DirNotFound');
if (autodiff_type == 1)
    rmpath(genpath(path_adi_2));
    addpath(genpath(path_adi_1));
elseif (autodiff_type == 2)
    rmpath(genpath(path_adi_1));
    addpath(genpath(path_adi_2));
end

addpath('/Users/sanderse/Dropbox/work/Programming/libs/');
addpath('RKproperties/');
addpath(genpath('/Users/sanderse/Dropbox/work/Programming/libs/chebfun-master/'));


%% select case file

% example case names (see case_files directory):
% Brusselator, Burgers, diffusion, Lorenz63, LotkaVolterra, polynomial,
% ProtheroRobinoson, reaction, simple, twobody, vanderpol

folder_cases = 'case_files';
case_name    = 'diffusion_filter';

folder_name  = [folder_cases '/' case_name '/'];
disp(['testcase: ' case_name]);

%% select RK method

% choose a RK method (see getRKmethod for options)
% implicit: GL1, GL2, GL3, RIIA1 (=CE), RIIA2, RIIA3, CHDIRK3, CHCONS3
% explicit: FE11, SSP22, RK44, Wray3
% dirk: SDIRK34
% esimplicit: CN22, LIIIA3
% esdirk: CN22, CHDIRK3

% RK_list = {'RK44'};
RK_list = {'GL1'};
% RK_list = {'CHC3','CHC5','GL2','GL3','RIIA2','RIIA3'};
RK_number = length(RK_list);


%% load case and settings
warning('off','all');
rmpath(genpath('case_files/'));
warning('on','all');
addpath(folder_name);
run([folder_name '/' case_name '.m']);

options.model.M  = M;           % dimension of the problem
options.model.u0 = u_start;     % initial condition

% number of simulations for time step refinement study (total number is
% Nsim * RK_number)
Nsim = length(N_list);

% error vector
err  = zeros(Nsim,RK_number);


%% loop over different RK methods
for RK_method = 1:RK_number
    
    
    %% get RK method
    
    options.RK_method.name = RK_list{RK_method};
    
    % get RK tableau and properties
    disp('*****************************************');
    [A_RK,b_RK,c_RK,r] = getRKmethod(options.RK_method.name);
    RK_order           = check_orderconditions(A_RK,b_RK,c_RK);
    disp(['RK method: ' RK_list{RK_method} ', order ' num2str(RK_order)]);
    
    s                        = length(b_RK);
    options.RK_method.stages = s;
    options.RK_method.A_RK   = A_RK;
    options.RK_method.b_RK   = b_RK;
    options.RK_method.c_RK   = c_RK;
    
    % check if explicit, DIRK, or fully implicit
    options.RK_method.type = check_RK_type(A_RK);
    disp(options.RK_method.type)
    
    
    %% initialize some variables
    
    % set counter for number of F evaluations (this includes Jacobian building)
    rhs_evals = 0;
    
    %
    options.time.t_start = t_start;
    options.time.t_end   = t_end;
    
    %
    options.settings.eps_nonlinear = eps_nonlinear;
    options.settings.it_max        = it_max;
    if (~exist('jacobian','var'))
        jacobian  = 1;
    end
    options.settings.jacobian      = jacobian;
    if (~exist('newton','var'))
        newton = 2;
    end
    options.settings.newton      = newton;
    
    options.settings.autodiff_type = autodiff_type;
    
    if (options.settings.newton == 0)       
        % determine Jacobian only once and its LU decomposition
        % note: this requires linear F and constant time step
        Fname = [options.model.function '_F'];
        options.settings.J =  Jacobian_single(Fname,u_start,t_start,options);
        % form iteration matrix, which is fixed during time steps
    end
    
    %% start time integration
    
    % loop over different time step values (for time step convergence studies)
    for j=1:Nsim
        
        % number of time steps
        N = N_list(j);
        
        % initial condition
        un = u_start;
        tn = t_start;
        
        % create vectors to store solution data
        u_plot  = zeros(M,N+1);
        t_plot  = zeros(N+1,1);
        nsolves = zeros(N+1,1);
        res     = zeros(N+1,1);
        u_plot(:,1) = un;
        t_plot(1)   = tn;
        
        % time step
        dt    = (t_end-t_start)/N;
        
        disp('===============');
        disp(['dt=' num2str(dt)]);
        
        tic
        
        if (options.settings.newton == 0)
            I_sM  = kron(speye(s),speye(M));
            Q     = (I_sM/dt-kron(A_RK,options.settings.J));
            % determine LU decomposition
            [L,U] = lu(Q);
            % store LU
            options.settings.Jacobian_L = L;
            options.settings.Jacobian_U = U;
        end
        
        % loop from t_start to t_end
        for i=2:N+1
            
            
            switch options.RK_method.type
                
                
                case 'dirk'
                    
                    [unew, nsolves(i), res(i)] = RKstep_dirk(un,tn,dt,options);
                    
                case 'esdirk'
                    
                    [unew, nsolves(i), res(i)] = RKstep_esdirk(un,tn,dt,options);
                    
                case 'implicit'
                    
                    [unew, nsolves(i), res(i)] = RKstep_implicit(un,tn,dt,options);
                    
                case 'esimplicit'
                    
                    [unew, nsolves(i), res(i)] = RKstep_esimplicit(un,tn,dt,options);
                    
                case 'explicit'
                    
                    unew = RKstep_explicit(un,tn,dt,options);
                    
                otherwise
                    error('wrong RK type');
                    
            end
            
            % new time level
            tnew = tn + dt;
            
            % update old solution and time level
            un   = unew;
            tn   = tnew;
            
            % store results for plotting purposes
            u_plot(1:M,i) = unew;
            t_plot(i)     = tnew;
        end
        
        toc
        
        if ( abs(t_end - tnew) > 1e-12)
            warning('end time not exactly reached');
        end
        
        % get exact solution (if available) at end time
        [u_exact,flag] = exact_solutions(t_end,options);
        
        if (flag==1)
            err(j,RK_method) = max(abs(unew(:) - u_exact(:)));
        end
        
        disp(['total number of rhs evaluations ' num2str(rhs_evals)])
        
    end
    
end

%% postprocessing

pp_file = [folder_name case_name '_pp.m'];
if (exist(pp_file,'file'))
    run(pp_file);
else
    disp('no postprocessing file found');
end

