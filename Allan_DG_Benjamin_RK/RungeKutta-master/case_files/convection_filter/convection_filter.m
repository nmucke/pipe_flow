%% PDE:
% du/dt = -c*du/dx
% discretized as
% du/dt = A*u = -C*u
% with periodic boundary conditions

% filtering:
% ubar(x) = int G(x,y) u(y) dy

% for no filter, use H0=dx, type = 'uni', width=0, skip=1

% DF: discrete first, filter next (new approach)
% FD: filter first, discretize next (classic approach)
approach = 'FD'; 

% reconstruction:
% 'inv': inverse of filter matrix to reconstruct
% 'RR': build sparse matrix with reconstruction rule
reconstruct = 'RR'; 

%% physical / numerical values

L     = 2*pi; % length of domain
% alpha = 0; % diffusion coefficient
c     = 1; % convection coefficient

M     = 20; % dimension of the model = number of interior grid points

dx    = L/M;


%% filter
% filtering width H; note full width is 2*H
% - see also filter_width.m for filter type
% often it makes sense to connect this to the grid width
H0 = L/4; % independent of grid
% H0 = dx/2; % depending on grid
filter_type = 'cos'; % 'cos' or 'uni'


%% discretization
options.model.function     = 'convection_filter';
options.model.constants(1) = c; % convection coefficient
options.model.constants(2) = H0; % filter width

options.model.filter_type  = filter_type; % filter width

x       = linspace(0,L,M+1)';  % grid points including boundary
x_in    = x(1:end-1); % periodic bc, leave out last point
% dx      = (x(end)-x(1))/M; % grid spacing, this equals diff(x)

options.model.x_in = x_in;


%% initial condition
t_start = 0;
t_end   = 0.5;

N_list  = 100;
k = 1;
u_start = sin(k*x_in-t_start); % initial condition before filtering

options.model.constants(3) = k; % frequency of initial condition

% note: dt = (t_end - t_start)/N_list

% settings for solving nonlinear equation
eps_nonlinear = 1e-12;
it_max  = 10;
jacobian = 2; % 1: FD, 2: AD
newton  = 1; % would be nice to make a special case in the code for linear problems, where the Jacobian only has to be determined once (newton=0)


%% quadrature
if (strcmp(approach,'DF'))
    % number of points used to construct quadrature rule
    % width = 2; % degree = 2*width+1
    % the effective width of the quadrature rule should be related to the
    % filter width
    % note: stencil = i-width:skip:i+width

    % for H0=L/4:
    width = floor(0.5*H0/dx); %M/5;
    skip  = 2; %width;

    % for H0=dx:
%     width = 1; %floor(H0/dx); %M/5;
%     skip  = 1; %width;
    
    % general:
%     width = floor(H0/dx); %M/5;
%     skip  = width;    
end



%% filter; note width = 2*H
[H, Hx] = filter_width(x_in,H0,filter_type);

% filter volumes
xf = [x_in-H x_in+H];

%% discrete filtering
% the filter or weight matrix maps from point values to 
% averaged values (both have same dimensions)

if (strcmp(approach,'DF'))

    addpath('/Users/sanderse/Dropbox/work/Programming/UQ/quadrature_rules');

    % weight matrix
    W     = zeros(M,M);

    % diffusion coefficient due to filter non-uniformity
    % alpha = zeros(M,1);

    % extend spatial coordinates for periodicity
    x_ext = [x_in-L; x_in; x_in+L];

    for i=1:M

        % stencil
        stencil = i-width:skip:i+width;
        % index in matrix
        indx    = mod(stencil,M);
        indx(indx==0) = M;

        % index in extended coordinates
        i_ext = i+M;
        indx2 = i_ext-width:skip:i_ext+width;
        % point set used:
        x0    = x_ext(indx2); 

        % integration domain
        domain = xf(i,:); %[x_in(i)-H(i),x_in(i)+H(i)];

        % row in weight matrix:
        W(i,indx) = QRgetWeights(length(x0),x0,domain);
        %

    end
    if (min(min(W))<0)
        warning('negative weights in W');
    end
    W = sparse(W);

elseif (strcmp(approach,'FD'))
    
    % additional term in PDE (from Taylor expansion)
    % diffusion:
    alpha    = (1/3).*H.*Hx;
    diagonal = ones(M,1);
    D        = spdiags(alpha,0,M,M)*(1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],M,M);
    D(1,end) = D(1,2);
    D(end,1) = D(end,end-1);
    
    
end

%% construct inverse of weighting matrix: 'reconstruction'
% the reconstruction matrix maps from averages to point values (both having
% same dimension)

if (strcmp(approach,'DF') && strcmp(reconstruct,'RR'))

    % reconstruction weight matrix
    R     = zeros(M,M);

    % extend spatial coordinates for periodicity
    x_ext = [x_in-L; x_in; x_in+L];
    x_bnd = (x_ext(1:end-1) + x_ext(2:end))/2;
    xf_ext = [x_in-L-H x_in-L+H; x_in-H x_in+H; x_in+L-H x_in+L+H];
    for i=1:M

        % stencil ('midpoints') of averages under consideration
%         width = 1;
%         skip  = 2;
%         stencil = i-width:skip:i+width;
        
        % index used in matrix
%         indx    = mod(stencil,M);
%         indx(indx==0) = M;

        % index in extended coordinates; used for getting boundary
        % coordinates
%         i_ext = i+M;
%         indx2 = i_ext-width-1:skip:i_ext+width;
%         Nbnd  = length(indx2);
%         domain = [x_bnd(indx2(1)),x_bnd(indx2(end))];

        % find all finite volume midpoints that lie within averaging volume
        % in extended coordinates
        
%         indx3 = i_ext-width:skip:i_ext+width;
        indx2 = find( (x_ext>xf(i,1) & x_ext<xf(i,2)) ==1);
        % take subset
%         indx3 = indx2(1):skip:indx2(end);
        indx3 = [indx2(1) indx2(ceil(end/2)) indx2(end)];
%         width = floor(sum(indx3)/2);
%         stencil = i-width:i+width;
        indx    = mod(indx3,M);
        indx(indx==0) = M;
        
        % row in weight matrix:
%         R(i,indx) = RRgetWeights([],x_in(i),Nbnd);
        R(i,indx) = RRgetWeights(xf_ext(indx3,:),x_in(i));
        %

    end
    if (min(min(R))<0)
        warning('negative weights in R');
    end
    R = sparse(R);

end

%% discretization matrices

% convection:
diagonal = ones(M,1);
C        = c*(1/(2*dx))*spdiags([-diagonal diagonal],[-1 1],M,M);
C(1,end) = -C(1,2);
C(end,1) = -C(end,end-1);

if (strcmp(approach,'DF'))

    % filtered convection matrix:
    if (strcmp(reconstruct,'RR'))
        C_filtered = W*C*R;
    elseif (strcmp(reconstruct,'inv'))
        C_filtered = W*C/W;
    end
    
    options.model.A = -C_filtered;

    % filtered initial condition
    u_start = W*u_start; 
    
elseif (strcmp(approach,'FD'))   

    options.model.A = -C + D;
    
    % non-uniform filtered solution:
    
    u_start = -(cos(k*(x_in+H)-t_start) - cos(k*(x_in-H)-t_start))./(2*k*H);
    
    % dt estimate
    dt_est = dx^2/max(abs(alpha))
    % dx estimate
    dx_est = 2*min(abs(alpha))
    
end

% plotting:
% imshow(W,'colormap',parula,'InitialMagnification','fit')
% or:
% imagesc(W)

%% optional: check what discretization the filtered matrix corresponds to
% use fixed width dx
% (better would be to build the vandermonde matrix for each grid point
% based on the filter widths)
if (strcmp(approach,'DF'))
    
    [~,diags] = spdiags(C_filtered);
    diag_core = diags(abs(diags)<M/2); % remove contributions 'far away', which are likely periodic BC
    mindiag   = min(diag_core);
    maxdiag   = max(diag_core);
    stencil   = mindiag:maxdiag; 
    
    n_stencil = length(stencil);
    V         = zeros(n_stencil+2,n_stencil);
    % build extended vandermonde to also get high order terms of Taylor
    % expansion (this doesn't work with the vander command)
    for i=1:n_stencil+2
        V(i,:) = 1./factorial(i-1) * (stencil*dx).^(i-1);
    end
    
    % u_eff contains the terms of the Taylor expansions that are involved,
    % e.g. u_eff(2) represents coeff*du/dx, where coeff can include grid
    % size effects (best to use dx=1)
    u_eff     = V*C_filtered(ceil(n_stencil/2),1:n_stencil)'
    
end

%% optional: check polynomial accuracy of R-W pair
if (strcmp(approach,'DF') && strcmp(reconstruct,'RR'))

    [~,diags] = spdiags(R*W);
    diag_core = diags(abs(diags)<M/2); % remove contributions 'far away', which are likely periodic BC
    mindiag   = min(diag_core);
    maxdiag   = max(diag_core);
    stencil   = mindiag:maxdiag; 
    n_stencil = length(stencil);

    order_poly = 2;
    coeff_test = rand(1,order_poly+1); % polynomial coefficients
    u_test     = polyval(coeff_test,x_in);
    u_int      = polyint(coeff_test); % note: u_avg is then effectively 1 degree higher
    x_ext      = [x; x(1)];
    x_bnd      = (x_ext(1:end-1) + x_ext(2:end))/2;
    u_avg      = polyval(u_int,x_bnd(2:end))-polyval(u_int,x_bnd(1:end-1))/dx; 

    % disregard boundary conditions, as the tested function is not periodic
    indx_interior = maxdiag+1:M-maxdiag-1;
    u_reconstruct = R*W*u_test;
    max(abs(u_reconstruct(indx_interior) - u_test(indx_interior)))

    u_weighted    = W*R*u_avg;
    max(abs(u_weighted(indx_interior) - u_avg(indx_interior)))
end

%% optional: plot reconstructed solution based on initial condition
if (strcmp(approach,'DF') && strcmp(reconstruct,'RR'))
    plot(x_in,u_start); % this is already the filtered initial condition
    hold on
    plot(x_in,R*u_start);
    legend('filtered IC','reconstructed IC')
elseif (strcmp(approach,'DF') && strcmp(reconstruct,'inv'))
    plot(x_in,u_start); % this is already the filtered initial condition
    hold on
    plot(x_in,W\u_start);
    legend('filtered IC','reconstructed IC')
end
    
