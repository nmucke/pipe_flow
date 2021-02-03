%% in this script we test kernel filtering a given signal on both periodic and compact domain
% we focus on the continuous case, i.e. use accurate integration methods,
% in order to investigate properties like conservation and normalization


%%
clear all
close all


%% inputs

% domain of signal for the compact/periodic domain case
x0 = -1;
x1 = 1;
% number of points (mainly for plotting purposes)
M = 101; % should take odd value

% the signal to be filtered (should be periodic on periodic domains)
u_start =@(x) sin(2*pi*x) + sin(8*pi*x);

% kernel:
kernel_type = 'box'; % options: spectral, gauss, box
% domain over which filter is applied - typically equal to x0, x1
% in case of 'gauss', the filtering domain will be taken -inf, inf;
% independent of what is specified below; for 'spectral', the filtering
% domain will be -1e2, 1e2
ksi0 = x0;
ksi1 = x1;
Mksi = M;

kernel_plot = 11; % number of kernels to plot (approximately)


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


%% set up domain and kernel functions
x_grid   = linspace(x0,x1,M)';
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
        % kernel for box filter on periodic domain
        G =@(r) (1/Delta).*(abs(r)<=Delta/2);
        % check normalization condition
        integral(@(r) G(r),ksi0,ksi1)
        % the normalization condition is not equal to 1 in the presence of boundaries
        % therefore we introduce the following normalization factor:
        % K = @(x) integral(@(ksi) G(x(i)-ksi),ksi0,ksi1);
        % this integral can be explicitly evaluated as
        K = @(x) 1 - (0.5 - (x-x0)/Delta).*(x<=x0+Delta/2) - ( 0.5 - (x1-x)/Delta).*(x>=(x1-Delta/2));
        % therefore, the expression for the effective kernel of the box filter on a bounded domain is:
        H = @(x,ksi) (1./K(x)).*(1/Delta).*(abs(x-ksi)<=(Delta/2 + eps));
        
        
    otherwise
        error('wrong filter type');
end

% we can now evaluate the following filter property:
% conservation
C = @(ksi) integral(@(x) H(x,ksi),x0,x1);

% and construct a new, symmetric kernel like Vreman:
J = @(x,ksi) 0.5*(H(x,ksi) + H(ksi,x)) + 1/(x1-x0) - 0.5*(arrayfun(C,x) + arrayfun(C,ksi));


%% get filtered signals on a periodic or unbounded domain
% (convolution-type filter)

% we want to compute the following integral:
% u_{filter}(x) = int u(ksi)*G(x-ksi) d ksi
% there are basically two approaches to compute this

%    approach 1, simply compute the integral as given:
%     u_{filter}(x) = int u(ksi)*G(x-ksi) d ksi

%    approach 2:
%    rewrite in terms of u; this requires an infinite or periodic domain and
%    symmetry of the kernel
%    the advantage is that it is easy to evaluate u outside [ksi0,ksi1] since
%    it is a periodic function, whereas it is not easy to extend G outside
%    [ksi0,ksi1]
%     u_{filter}(x) = int u(x-ksi)*G(ksi) rd ksi

% the integrals over ksi are accurately approximated using integral(); note
% that we are not using ksi_grid here

u_filter1 = zeros(M,1);
u_filter2 = zeros(M,1);

for i=1:length(x_grid)
    
    
    switch kernel_type
        case 'gauss'
            % approach 1:
            u_filter1(i) = integral(@(ksi) u_start(ksi).*G(x_grid(i)-ksi),-inf,inf);
            % approach 2:
            u_filter2(i) = integral(@(ksi) u_start(x_grid(i) - ksi).*G(ksi),-inf,inf);
        case 'spectral'
            % approach 1
            u_filter1(i) = integral(@(ksi) u_start(ksi).*G(x_grid(i)-ksi),-1e2,1e2);
            % approach 2
            u_filter2(i) = integral(@(ksi) u_start(x_grid(i) - ksi).*G(ksi),-1e2,1e2);
        case 'box'           
            % approach 1
            u_filter1(i) = integral(@(ksi) u_start(ksi).*G(x_grid(i)-ksi),ksi0-Lksi,ksi1+Lksi);
            % approach 2
            u_filter2(i) = integral(@(ksi) u_start(x_grid(i) - ksi).*G(ksi),ksi0-Lksi,ksi1+Lksi);
            
    end
    
    
end


figure(1)
r_plot = x_grid;
plot(x_grid,u_start(x_grid));
hold on
plot(r_plot,G(r_plot));
plot(x_grid,u_filter1,'s-');
plot(x_grid,u_filter2,'x-');
legend('original','kernel','filtered, approach 1','filtered, approach 2');
        
title('periodic domain filtering')


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


%% plots
figure(2)
plot(x_grid,u_start(x_grid),'-');
hold on
plot(x_grid,u_bounded_filter1,'s-');
plot(x_grid,u_bounded_filter2,'s-');
plot(x_grid,u_bounded_filter3,'s-');
legend('original','filtered, approach 1','filtered, approach 2','filtered, approach 3');
title('compact domain filtering')

% plot effective kernels H and the new conservative kernel J
figure(3)
skip = max(floor((M-1)/kernel_plot),1);
for i=1:skip:M

    p1 = plot(ksi_grid,H(x_grid(i),ksi_grid),'-');
    hold on
    plot(ksi_grid,J(x_grid(i),ksi_grid),'--','Color',p1.Color);
end
legend('H(x,\xi)','J(x,\xi)');
% title('effective kernel J(x,\xi)')

% plot normalization factor K(x) and compare with Gnorm (should be the
% same)
% note: differences between G_norm2 and K can appear due to numerical integration
% errors
figure(4)
plot(x_grid,G_norm2)
hold on
plot(x_grid,K(x_grid),'s--')
title('normalization factor K(x)');

%% plot conservation error
figure(5)
p1 = plot(ksi_grid,Cons-1,'-');
hold on
plot(ksi_grid,Cons_J-1,'--','Color',p1.Color);
title('conservation error');
legend('H(x,\xi)','J(x,\xi)');
grid
xlabel('\xi')