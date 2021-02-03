%% in this script we test kernel filtering a given signal using discrete filters

clear all
close all
format long
format compact

addpath('/Users/sanderse/Dropbox/work/Programming/UQ/quadrature_rules');
addpath('/Users/sanderse/Dropbox/work/Programming/libs/chebfun-master/');


fontsize = 16;
fontname = 'Helvetica';
linewidth = 2;
set(0,'defaultlinelinewidth',linewidth)

%% inputs

% domain of signal for the compact/periodic domain case
x0 = 0;
x1 = 1;
% number of points
Mfine = 201; % for plotting fine solution, e.g. original signal
M_list = [11 21 41]; % 81 161]; % should take odd value


% the signal to be filtered (should be periodic on periodic domains)
u_start =@(x) sin(2*pi*x).^2 + sin(8*pi*x);

% kernel:
kernel_type = 'box'; % options: spectral, gauss, box
% domain over which filter is applied - typically equal to x0, x1
% in case of 'gauss', the filtering domain will be taken -inf, inf;
% independent of what is specified below; for 'spectral', the filtering
% domain will be -1e2, 1e2
ksi0 = x0;
ksi1 = x1;
Mksi_list = M_list;

kernel_plot = 11; % number of kernels to plot (approximately)

% quadrature rule parameters
QRtype_local = 'trapz'; % options: 'trapz','poly','poly_kernel'
QRtype_global = 'trapz'; % options: 'trapz','poly'

% kernel parameters:
switch kernel_type
    case 'gauss'
        gamma = 6;
        Delta = 0.2;
    case 'spectral'
        % specify cut off frequency:
        k_c = 4*pi;
        % base Delta on cut off:
        Delta = pi/k_c;
    case 'box'
        % note that the transfer function is such that G(k*Delta)=0 for
        % k*Delta = 2*pi, 4*pi, 6*pi, 8*pi, etc.
        % so for e.g. Delta = 1/2, we get that k=4*pi, 8*pi, 12*pi etc are
        % completely damped
        Delta = 0.2;
end


%%
eps = 1e-16; % machine precision value used to avoid issues at r=0
Nsims = length(M_list);
quad_error = zeros(Nsims,1);

for q=1:length(M_list)
    
    M = M_list(q);
    Mksi = Mksi_list(q);
    
    %% set up domain and kernel functions
    x_grid   = linspace(x0,x1,M)';
    x_fine   = linspace(x0,x1,Mfine)';
    Lx       = x1-x0;
    dx       = Lx/M;
    ksi_grid = linspace(ksi0,ksi1,Mksi)';
    Lksi     = ksi1-ksi0;
    dksi     = Lksi/Mksi;
    
    
    % we want to compute the following integral (a kernel filter):
    % u_{filter} = int u(ksi)*H(x,ksi) d ksi
    
    % on an unbounded or periodic domain, with uniform filter size, the integral is of convolution type:
    % u_{filter} = int u(ksi)*G(x-ksi) d ksi
    % in this case the kernel only depends on r = x - ksi, so we have
    % G(x-ksi) = H(x,ksi)
    
    % kernel, assumed to decay fast enough to be used in periodic integration
    % on bounded domains the function is truncated and scaled such to integrate
    % to 1
    switch kernel_type
        case 'gauss'
            G =@(r) sqrt( (gamma/(pi* (Delta^2)))) * exp(- gamma*(r.^2) / Delta^2);
            % check if G integrates to one (normalization condition)
            integral(@(r) G(r),-inf,inf)
            
            % we also introduce the more generic filter kernel, which does not
            % assume convolution property:
            H = @(x,ksi) sqrt( (gamma/(pi* (Delta^2)))) * exp(- gamma*(x-ksi).^2 / Delta^2);
            
        case 'spectral'
            G =@(r) (sin(pi*(r+eps)/Delta)./(pi*(r+eps)));
            % check if G integrates to one (normalization condition)
            % note: infinity gives issues in accuracy
            integral(@(r) G(r),-1e2,1e2)
            
            % we also introduce the more generic filter kernel, which does not
            % assume convolution property:
            H = @(x,ksi) sin(pi*(x-ksi)/Delta)./(pi*(x-ksi));
            
        case 'box'
            G =@(r) (1/Delta).*(abs(r)<=Delta/2);
            % check normalization condition
            integral(@(r) G(r),ksi0,ksi1)
            % the normalization condition is not equal to 1 in the presence of boundaries
            % therefore we introduce the following normalization factor:
            % K = @(x) integral(@(ksi) G(x(i)-ksi),ksi0,ksi1);
            % this integral can be explicitly evaluated as
            K = @(x) 1 - (0.5 - (x-x0)/Delta).*(x<=x0+Delta/2) - ( 0.5 - (x1-x)/Delta).*(x>=(x1-Delta/2));
            % find the integration bounds where G is nonzero,
            % these are such that K = D_right - D_left
            % these will be used to set up the quadrature rules
            D_left  = @(x) - 0.5 + (0.5 - (x-x0)/Delta).*(x<=x0+Delta/2);
            D_right = @(x) 0.5 - ( 0.5 - (x1-x)/Delta).*(x>=(x1-Delta/2));
            % expression for effective kernel
            H = @(x,ksi) (1./K(x)).*(1/Delta).*(abs(x-ksi)<=(Delta/2 + eps));
            
        otherwise
            error('wrong filter type');
    end
    
    % we can now evaluate the following filter property:
    % conservation
    C = @(ksi) integral(@(x) H(x,ksi),x0,x1);
    
    % and construct a new, symmetric kernel like Vreman:
    J = @(x,ksi) 0.5*(H(x,ksi) + H(ksi,x)) + 1/(x1-x0) - 0.5*(arrayfun(C,x) + arrayfun(C,ksi));
    
    
    %% get filtered signal on a non-periodic, bounded domain
    % note: this only makes sense for the box filter
    
    % the convolution approach outlined above does not apply in this case, as
    % the kernel G(x-ksi) depends on the spatial location
    % we therefore need to use the more generic form H(x,ksi):
    
    % u_{filter}(x) = int u(ksi)*H(x,ksi) d ksi
    
    u_bounded_filter1 = zeros(M,1);
    u_bounded_filter2 = zeros(M,1);
    u_bounded_filter3 = zeros(M,1);
    
    G_norm2 = zeros(M,1);
    
    for i=1:length(x_grid)
        
        % simple integration
        G_bounded = G(x_grid(i) - ksi_grid);
        G_norm1   = sum(G_bounded.*dx);
        u_bounded_filter1(i) = sum(u_start(ksi_grid).*G_bounded.*dx)/G_norm1;
        
        % more accurate integration:
        % note that G_norm2 should equal K(x) (this is tested below)
        G_norm2(i) = integral(@(ksi) G(x_grid(i)-ksi),ksi0,ksi1);
        u_bounded_filter2(i) = integral(@(ksi) u_start(ksi).*G(x_grid(i)-ksi),ksi0,ksi1)/ G_norm2(i);
        
        % alternatively, we can directly integrate H(x,ksi), this should give
        % an equivalent result
        u_bounded_filter3(i) = integral(@(ksi) u_start(ksi).*H(x_grid(i),ksi),ksi0,ksi1);
    end
    
    % conservation tests
    Cons   = zeros(M,1);
    Cons_J = zeros(M,1);
    %
    for i=1:length(ksi_grid)
        
        Cons(i)   =  C(ksi_grid(i)); %integral(@(x) H_box(x,ksi_grid(i)),x0,x1);
        % this is expensive, since J itself also includes an integral
        %     Cons_J(i) =  integral(@(x) J(x,ksi_grid(i)),x0,x1);
        
    end
    
    %     figure(2)
    %     plot(x_fine,u_start(x_fine),'-');
    %     hold on
    % plot(x_grid,u_bounded_filter1,'s-');
    % plot(x_grid,u_bounded_filter2,'s-');
    %     plot(x_grid,u_bounded_filter3,'s-');
    %     title('compact domain filtering')
    
    %% we now proceed with discretely approximating the filtering operation by using a quadrature rule
    % first, we set up a quadrature rule to simply integrate a function over the domain
    % without weights that are specific to the kernel,
    % in other words the filtered signal will be u_{filtered} = (W.*G)*u, with W
    % the weights of the quadrature rule
    
    % set-up weight matrix
    W_trapz = zeros(M,M);
    W_poly  = zeros(M,M);
    W       = zeros(M,M);

    
    % ratio of filter length to mesh width, typically larger than 1
    fg_ratio = Delta / dx;
    if (fg_ratio < 2)
        error('filter-grid ratio < 2');
    end
    % width of integration rule; number of points is 2*width+1
    % nwidth   = 2; % degree = 2*width
    
    % allocate kernel
    H_d  = zeros(M,Mksi);
    
    for i=1:M
        
        % integration domain
        x_left  = max(x_grid(i) - Delta/2,x0);
        x_right = min(x_grid(i) + Delta/2,x1);
        
        domain = [x_left,x_right];
        dx_domain = diff(domain);
        
        % stencil
        %     i_left  = max(i-nwidth,1);
        %     di_left = i_left - (i-nwidth);
        %     i_right = min(i+nwidth,M);
        %     di_right = (i+nwidth) - i_right;
        %     stencil  = i_left-di_right:i_right+di_left;
        
        % assume nodes are chosen that lie within the integration domain (this
        % is not necessary but is typically more accurate)
        stencil = find(x_grid>=(x_left-eps) & x_grid<=(x_right+eps));
        x_nodes = x_grid(stencil);
        Mi      = length(x_nodes);
        
        switch QRtype_local
            
            case 'trapz'
                % local integration with trapezoidal rule
                % weighting with kernel performed afterwards
                % row in weight matrix:
                % note that QRgetWeights returns the average so this needs to be
                % scaled:
                W_trapz(i,stencil) = dx_domain*QRgetWeights(Mi,x_nodes,domain,[],[],[],[],'trapz');
                
                % kernel
                H_d(i,:) = H(x_grid(i),ksi_grid);

                % effective weight matrix
                W(i,:) = W_trapz(i,:).*H_d(i,:);
                
            case 'poly'
                
                % local integration with high order polynomial
                W_poly(i,stencil) = dx_domain*QRgetWeights(Mi,x_nodes,domain);
                
                % kernel
                H_d(i,:) = H(x_grid(i),ksi_grid);

                % effective weight matrix
                W(i,:) = W_poly(i,:).*H_d(i,:);
                
            case 'poly_kernel'
                
                % do weighting with the kernel in constructing W
                distribution.type = 'unif';
                distribution.x_left = x_left;
                distribution.x_left = x_right;
                % note that here we are not multiplying by dx_domain as we build
                % directly the approximation to the filtered solution
                W(i,stencil) = QRgetWeights(Mi,x_nodes,domain,distribution)';
        end
    end
    
    ubar_quad = W*u_start(x_grid);
    
    % check normalization: right multiplication of W
    norm_error(q,1) = max(abs(ones(M,1) - W*ones(M,1)));
    disp(['normalization error: ' num2str(norm_error(q))]);
    
    % we also need a quadrature rule to integrate over the entire domain,
    % in order to investigate conservation properties
    % weight matrix
    %     A_d   = zeros(1,M);
    % stencil and point selection
    stencil = 1:M;
    x_nodes = x_grid(stencil);
    % integration domain
    domain = [x0,x1];
    
    switch QRtype_global
        case 'trapz'
            % a single row in the weight matrix can be set up by calling the QR routine, e.g.
            A_d = QRgetWeights(M,x_nodes,domain,[],[],[],[],'trapz')';
            % where x_nodes are the nodes to be used in the quadrature rule, and domain
            % is the domain of integration.
        case 'poly'
            % this is likely to be inaccurate due to negative weights with high order polynomial
            % quadrature on equidistant nodes
            A_d  = QRgetWeights(M,x_nodes,domain);
            
    end
    
    % check conservation properties: left multiplication of W
    cons_error_sol(q,1) =  max(abs((A_d * W  - A_d)*u_start(x_grid)));
    cons_error_vnorm(q,1) =  max(abs(A_d * W  - A_d));
    disp(['conservation error solution: ' num2str(cons_error_sol(q))]);
    disp(['conservation error vector norm: ' num2str(cons_error_vnorm(q))]);
    
    % make conservative:
    % require that A_d * Wnew * u_start = A_d * u_start
    % (for any u!)
    e    = ones(M,1);
    T    = (1/sum(A_d)) * kron(e,A_d*(speye(M) - W));
    Wnew = W + T;
    ubar_cons = Wnew*u_start(x_grid);
    
    % check conservation properties: left multiplication of Wnew
    cons_error_sol_cons(q,1) =  max(abs((A_d * Wnew  - A_d)*u_start(x_grid)));
    cons_error_vnorm_cons(q,1) =  max(abs(A_d * Wnew  - A_d));
    disp(['conservation error solution: ' num2str(cons_error_sol_cons(q))]);
    disp(['conservation error vector norm: ' num2str(cons_error_vnorm_cons(q))]);
    
    figure(1)
    plot(x_grid,ubar_quad,'d-');
    hold on
    plot(x_grid,ubar_cons,'o--');
    %     legend('original','filtered','non-uniform filtered','non-uniform quadrature');
    
    quad_error(q,1) = max(abs(ubar_quad - u_bounded_filter3));
    quad_error_cons(q,1) = max(abs(ubar_cons - u_bounded_filter3));
    
end

%%
figure(1);
% set(sp,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
% set(sp,'defaultTextFontName', fontname)
set(gca,'FontSize',fontsize)

hold on
plot(x_grid,u_bounded_filter3,'-');
plot(x_fine,u_start(x_fine),'-');
xlabel('x')
ylabel('u')
grid on
legend([string(M_list),'reference','unfiltered'])

figure(2)
loglog(M_list,quad_error,'s-');
hold on
loglog(M_list,norm_error,'s-');
loglog(M_list,cons_error_sol,'s-');
loglog(M_list,cons_error_vnorm,'s-');
loglog(M_list,quad_error_cons,'s--');
loglog(M_list,cons_error_sol_cons,'s--');
loglog(M_list,cons_error_vnorm_cons,'s--');
xlabel('M');
ylabel('error');
legend('quadrature error','normalization error','conservation error solution',...
    'conservation error vector norm','quadrature error conservative',...
    'conservation error solution cons','conservation error vector norm cons');

xlim([10 200])
grid on
set(gca,'FontSize',fontsize)


%% plot conservation error
figure(5)
p1 = plot(ksi_grid,Cons-1,'-');
hold on
plot(ksi_grid,Cons_J-1,'--','Color',p1.Color);
title('conservation error');
legend('H(x,\xi)','J(x,\xi)');
grid
xlabel('\xi')