%% PDE:
% du/dt = -cx*du/dx - cy*du/dy + alpha*(d2u/dx^2 + d2u/dy^2)
% discretized as
% du/dt = A*u = -C*u + D*u + ybc

addpath('/Users/sanderse/Dropbox/work/Programming/UQ/quadrature_rules');

%% physical values
Lx    = 1; % length of domain
Ly    = 1; % length of domain
alpha = 0.01; % diffusion coefficient
cx    = 1; % convection coefficient
cy    = 0.25; % convection coefficient

%% averaging settings
average = 1;
average_type = 1; % 0: composite midpoint; 1: global polynomial; this is also used in case average=0, for postprocessing
My_new  = 3; % number of averaging volumes in y-direction used in averaging and reconstruction; only used when average=1

%% discretization
Mx     = 30; % dimension of the model = number of interior grid points
My     = 30; % dimension of the model = number of interior grid points
options.model.function     = 'convectiondiffusion2D';
options.model.constants(1) = alpha; % diffusion coefficient
% options.model.constants(2) = L; % domain length
options.model.constants(3) = cx; % convection coefficient
options.model.constants(4) = cy; % convection coefficient
options.model.bc.left.type  = 'dir';
options.model.bc.right.type = 'neu';
options.model.bc.low.type  = 'dir';
options.model.bc.up.type   = 'neu';

x       = linspace(0,Lx,Mx+1)';  % grid lines including boundary
y       = linspace(0,Ly,My+1)';  % grid lines including boundary
x_in    = (x(1:end-1) + x(2:end))/2; % grid points
y_in    = (y(1:end-1) + y(2:end))/2; % grid points

options.model.bc.left.val  = -2+2*cos(pi*y_in)-2*y_in.*(y_in.^2-3); %sin(2*pi*y_in).^2-5*y_in.*(y_in-1).^3; %-6*y_in.*(y_in-1);   %sin(2*pi*y_in).^2;
options.model.bc.right.val = 0;
options.model.bc.low.val  = 0;
options.model.bc.up.val   = 0;


dx      = (x(end)-x(1))/Mx; % grid spacing, this equals diff(x)
dy      = (y(end)-y(1))/My; % grid spacing, this equals diff(x)

options.model.x_in= x_in;
options.model.y_in= y_in;

%% postprocessing
options.pp.movie = 1;
options.pp.moviename = 'cospoly';
options.pp.movierate = 16; % frames per second


%% make BC vectors if necessary
if (length(options.model.bc.left.val) == 1)
    options.model.bc.left.val  = options.model.bc.left.val * ones(My,1);
end
if (length(options.model.bc.right.val) == 1)
    options.model.bc.right.val  = options.model.bc.right.val * ones(My,1);
end
if (length(options.model.bc.low.val) == 1)
    options.model.bc.low.val  = options.model.bc.low.val * ones(Mx,1);
end
if (length(options.model.bc.up.val) == 1)
    options.model.bc.up.val  = options.model.bc.up.val * ones(Mx,1);
end

%%
convection_diffusion_operators;

%% construct averaging and reconstruction operators
% the operators are also constructed when averaging is not used, in order
% to apply them in a postprocessing step when necessary

% % averaging over total domain
% diagonal   = ones(Mx,1);
% Av1D_total = (1/My)*ones(1,My);
% Ave_total  = kron(Av1D_total,spdiags(diagonal,0,Mx,Mx)); % ones could be replaced by a suitable quadrature rule

%% averaging

if (average == 1)
    
    Msub = My/My_new; %number of points per averaging volume
    
    switch average_type
        
        case 0
            % composite midpoint
            QR1D = 1/Msub*ones(1,Msub);
            
        case 1
            % N-degree polynomial QR
            QR1D = QRgetWeights(Msub,y_in(1:Msub),[0,1/My_new])';
            
    end
    
    % extend the QR to all averaging volumes
    Av1D = kron(speye(My_new),QR1D);
    
    % extend to 2D
    diagonal = ones(Mx,1);
    Ave  = kron(Av1D,spdiags(diagonal,0,Mx,Mx));
    
elseif (average == 0)
    % averaging operator is still constructed to average the result
    % (postprocessing only)
    
    switch average_type
        
        case 0
            % composite midpoint
            Av1D = 1/My*ones(1,My);
            
        case 1
            % N-degree polynomial QR
            Av1D = QRgetWeights(My,y_in,[0,1])';
            
    end
    
    % extend to 2D
    diagonal = ones(Mx,1);
    Ave  = kron(Av1D,spdiags(diagonal,0,Mx,Mx));
    
end



%% reconstruction

if (average == 1)
    % assume polynomial profile, e.g. u = a0*y^4 + a1*y^3 + a2*y^2 + a3*y + a4;
    % the degree is given by r = My_new+2
    % form Vandermonde, with My_new averaging conditions and two boundary
    % conditions: V=[1 1 1 1 1;0 0 0 0 1; etc.]; a = A\[0 0 1];
    r = My_new+2;
    
    % equidistant choice of averaging volumes:
    y_bnd = linspace(0,1,My_new+1)'; %[0;1/3;2/3;1];
    
    % polynomial basis 
    pol   = ones(1,r);
    p_int = polyint(pol);
    p_der = polyder(pol);
    
    % Vandermonde matrix:
    V = zeros(r,r);
    switch options.model.bc.up.type
        case 'dir'
            V(r-1,:) = ones(1,r); % Dirichlet condition upper
        case 'neu'
            arr  = 1:r-1;   
            for j=2:r
                arr_test = arr;
                arr_test(j-1) = [];                
                p_test = p_der;
                p_test(arr_test) = 0;
                V(r-1,j) = polyval(p_test(j-1),Ly); % Neumann condition upper
            end
    end
    switch options.model.bc.low.type
        case 'dir'
            V(r,end) = 1; % Dirichlet condition lower
        case 'neu'
            for j=2:r
                V(r,j) = polyval(p_der(j-1),0); % Neumann condition lower
            end   
    end

    arr  = 1:r;
    for i=1:r-2
        for j=1:r
            arr_test = arr;
            arr_test(j) = [];
            p_test   = p_int;
            p_test(arr_test) = 0;
            V(i,j) = (polyval(p_test,y_bnd(i+1)) - polyval(p_test,y_bnd(i)))/(y_bnd(i+1)-y_bnd(i));
        end
    end
    %     mu = [1;1;1;0;0];
    %     a  = V\mu;
    W = inv(V);
    
    % evaluate the reconstruction polynomial at the grid points
    phi_basis = eye(r);
    V_eval  = zeros(My,r);
    
    for i=1:r
        V_eval(:,i)  = polyval(phi_basis(i,:),y_in);
    end
    
    Rec1D = V_eval*W(:,1:r-2); % here we use that the BC are 0; if not zero, then an additional term will result that enters into a BC-vector
    diagonal = ones(Mx,1);
    Rec   = kron(Rec1D,spdiags(diagonal,0,Mx,Mx)); % ones could be replaced by a suitable quadrature rule
    
    
    % check if A*R = I
    % 1D:
    if (max(max(abs(Av1D*Rec1D-speye(My_new))))>1e-12)
        % full 2D: max(max(abs(Ave*Rec - speye(My_new*Mx))))>1e-12)
        warning('averaging and reconstruction operators not consistent');
    end
    
    %
    Cx_orig = Cx;
    Cy_orig = Cy;
    Dx_orig = Dx;
    Dy_orig = Dy;
    
    % construct new operators
    Cx = Ave*Cx*Rec;
    Cy = Ave*Cy*Rec;
    Dx = Ave*Dx*Rec;
    Dy = Ave*Dy*Rec;
    
    bCx = Ave*bCx;
    bCy = Ave*bCy;
    bDx = Ave*bDx;
    bDy = Ave*bDy;

    % error estimates:
    max(max(abs(Rec*Ave-speye(Mx*My))))
    
    max(max(abs(Cx*Ave-Ave*Cx_orig)))
    max(max(abs(Cy*Ave-Ave*Cy_orig)))
    max(max(abs(Dx*Ave-Ave*Dx_orig)))
    max(max(abs(Dy*Ave-Ave*Dy_orig)))
end

%% assemble: du/dt = -C*u + D*u + bc

options.model.A    = - Cx - Cy + Dx + Dy;
options.model.bc.r = - bCx - bCy + bDx + bDy;


%% initial condition
% u_start = sin(2*pi*x_in); % initial condition
% u_start = kron(ones(My,1),u_start);
u_start = -2+2*cos(pi*y_in)-2*y_in.*(y_in.^2-3); %sin(2*pi*y_in).^2-5*y_in.*(y_in-1).^3; %-6*y_in.*(y_in-1); %sin(2*pi*y_in).^2;  %1+sin(2*pi*y_in).^2; % initial condition
u_start = kron(u_start,ones(Mx,1));

if (average == 1)
    u_start = Ave*u_start;
end

M  = length(u_start);

%% time integration settings
t_start = 0;
t_end   = 1;

N_list  = 400;

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 1; % would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)