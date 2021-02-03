%% diffusion:

% in x-direction:
diagonal = ones(Mx,1);
Dx_1D     = (1/(dx^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],Mx,Mx);
bDx_left  = zeros(Mx,1);
bDx_right = zeros(Mx,1);

% boundary conditions, such that
% u_xx = Dx*u + bDx

switch options.model.bc.left.type
    case 'dir'
        Dx_1D(1,1)     = -3/(dx^2);
        bDx_left(1) = 2/(dx^2);
        
    case 'neu'
        Dx_1D(1,1)     = -1/(dx^2);
        Dx_1D(1,2)     = 1/(dx^2);
        bDx_left(1) = 1/dx;
        
    case 'per'
        Dx_1D(1,end)   = 1/(dx^2);
        
    otherwise
        error('wrong BC specification');
end

switch options.model.bc.right.type
    case 'dir'
        Dx_1D(end,end)    = -3/(dx^2);
        bDx_right(end) = 2/(dx^2);
        
    case 'neu'
        Dx_1D(end,end-1)  = 1/(dx^2);
        Dx_1D(end,end)    = -1/(dx^2);
        bDx_right(end) = 1/dx;
        
    case 'per'
        Dx_1D(end,1) = 1/(dx^2);
        
    otherwise
        error('wrong BC specification');
end

Dx_1D  = alpha*Dx_1D;
bDx_left  = alpha*bDx_left;
bDx_right = alpha*bDx_right;

% extend in y-direction:
Dx       = kron(speye(My),Dx_1D);
bDx      = kron(options.model.bc.left.val,bDx_left) + ...
           kron(options.model.bc.right.val,bDx_right);


% in y-direction:
diagonal = ones(My,1);
Dy_1D    = (1/(dy^2))*spdiags([diagonal -2*diagonal diagonal],[-1 0 1],My,My);
bDy_low  = zeros(My,1);
bDy_up   = zeros(My,1);

% boundary conditions, such that
% u_yy = Dy*u + bDy

switch options.model.bc.low.type
    case 'dir'
        Dy_1D(1,1)    = -3/(dy^2);
        bDy_low(1) = 2/(dy^2);
        
    case 'neu'
        Dy_1D(1,1)    = -1/(dy^2);
        Dy_1D(1,2)    = 1/(dy^2);
        bDy_low(1) = 1/dy;
        
    case 'per'
        Dy_1D(1,end)  = 1/(dy^2);
        
    otherwise
        error('wrong BC specification');
end

switch options.model.bc.up.type
    case 'dir'
        Dy_1D(end,end) = -3/(dy^2);
        bDy_up(end) = 2/(dy^2);
        
    case 'neu'
        Dy_1D(end,end-1) = 1/(dy^2);
        Dy_1D(end,end)   = -1/(dy^2);
        bDy_up(end)   = 1/dy;
        
    case 'per'
        Dy_1D(end,1) = 1/(dy^2);
        
    otherwise
        error('wrong BC specification');
end

Dy_1D  = alpha*Dy_1D;
bDy_low = alpha*bDy_low;
bDy_up  = alpha*bDy_up;

% extend in x-direction:
Dy       = kron(Dy_1D,speye(Mx));
bDy      = kron(bDy_low,options.model.bc.low.val) + ...
           kron(bDy_up,options.model.bc.up.val);

%% convection

diagonal = ones(Mx,1);
Cx_1D    = (0.5/dx)*spdiags([-diagonal diagonal],[-1 1],Mx,Mx);
bCx_left  = zeros(Mx,1);
bCx_right = zeros(Mx,1);


switch options.model.bc.left.type
    case 'dir'
        Cx_1D(1,1)     = 0.5/dx;
        Cx_1D(1,2)     = 0.5/dx;
        bCx_left(1) = -1/dx;
        
    case 'neu'
        Cx_1D(1,1)     = -0.5/dx;
        Cx_1D(1,2)     = 0.5/dx;
        bCx_left(1) = 0.5;
        
    case 'per'
        Cx_1D(1,end) = -0.5/dx;
        
    otherwise
        error('wrong BC specification');
end

switch options.model.bc.right.type
    case 'dir'
        Cx_1D(end,end)    = -0.5/dx;
        Cx_1D(end,end-1)  = -0.5/dx;
        bCx_right(end) = 1/dx;
        
    case 'neu'
        Cx_1D(end,end-1)  = -0.5/dx;
        Cx_1D(end,end)    = 0.5/dx;
        bCx_right(end) = 0.5;
        
    case 'per'
        Cx_1D(end,1) = 0.5/dx;
        
    otherwise
        error('wrong BC specification');
end

Cx_1D = cx*Cx_1D;
bCx_left = cx*bCx_left;
bCx_right = cx*bCx_right;

% extend in y-direction:
Cx  = kron(speye(My),Cx_1D);
bCx = kron(options.model.bc.left.val,bCx_left) + ...
      kron(options.model.bc.right.val,bCx_right);


% in y-direction:
diagonal = ones(My,1);
Cy_1D    = (0.5/dy)*spdiags([-diagonal diagonal],[-1 1],My,My);
bCy_low  = zeros(My,1);
bCy_up   = zeros(My,1);

switch options.model.bc.low.type
    case 'dir'
        Cy_1D(1,1) = 0.5/dy;
        Cy_1D(1,2) = 0.5/dy;
        bCy_low(1)  = -1/dy;
        
    case 'neu'
        Cy_1D(1,1) = -0.5/dy;
        Cy_1D(1,2) = 0.5/dy;
        bCy_low(1)  = 0.5;
        
    case 'per'
        Cy_1D(1,end) = -0.5/dy;
        
    otherwise
        error('wrong BC specification');
end

switch options.model.bc.up.type
    case 'dir'
        Cy_1D(end,end)   = -0.5/dy;
        Cy_1D(end,end-1) = -0.5/dy;
        bCy_up(end)   = 1/dy;
        
    case 'neu'
        Cy_1D(end,end-1)  = -0.5/dy;
        Cy_1D(end,end)    = 0.5/dy;
        bCy_up(end)    = 0.5;
        
    case 'per'
        Cy_1D(end,1) = 0.5/dy;
        
    otherwise
        error('wrong BC specification');
end

Cy_1D = cy*Cy_1D;
bCy_low = cy*bCy_low;
bCy_up  = cy*bCy_up;

% extend in y-direction:
Cy  = kron(Cy_1D,speye(Mx));
bCy = kron(bCy_low,options.model.bc.low.val) + ...
    kron(bCy_up,options.model.bc.up.val);
