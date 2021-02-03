% Trefethen results as in SIAM Review article
clearvars
close all


d = [-0.5,1];
testcase = 'jump1';

Mlist = [1; 2; 4; 8];

linestyle = {'-','-.','--',':','-','-.','--',':'};
fontsize = 14;
fontname = 'Helvetica'; 
linewidth = 3;

save_pdf = 0;

legendstr = {};

addpath(genpath('/Users/sanderse/Dropbox/work/Programming/libs/'))
for k=1:length(Mlist)
    
    M = Mlist(k); % subdivision of domain into M intervals
    
    switch testcase
        case 'Gaussian'
            lambda = 50;
            gg = @(x) exp(-lambda*x.^2);
        case 'Runge'
            gg = @(x) 1./(1+16*x.^2);
        case 'polynomial'
            gg = @(x) x.^20;
        case 'abs'
            gg = @(x) abs(x).^3;
        case 'dimple'
            gg = @(x) exp(-1./(x.^2));
        case 'jump0'
            gg = @(x) sin(x).*(x<0) + cos(x).*(x>0);
%             gg_ex = -(cos(0)-cos(-0.5)) + (sin(1)-sin(0));            
        case 'jump1'
            gg = @(x) sin(x).*(x<0) + -0.5*sin(x).*(x>0);
        otherwise
            error('testcase unknown');
    end
    % domain
    dx = diff(d)/M;
    % x = nonuniform_grid(dx,d(1),d(2),1.1);
    sx = 1;
    x=nonuniform_grid2(d(1),d(2),sx,M);
    % x = d(1):dx:d(2);
    % M = length(x)-1;
    
    % exact value from integrating the chebfun
    I = sum(chebfun(gg,d,'splitting','on'));
%     I =-(cos(0)-cos(-1)) + (sin(1)-sin(0));
    
    errcc = [];
    errgauss = [];
    nn = 1:1:60/M;
    for n = nn
        Icc = 0;
        Igauss = 0;
        
        for m=1:M
            dsub = [x(m),x(m+1)];
            
            % Clenshaw-Curtis integration = chebfun of degree n+1
            Icc = Icc + sum(chebfun(gg,dsub,n+1));
            
            % Gauss-Legendre points
            [s,w] = legpts(n+1,dsub);
            Igauss = Igauss + w*gg(s);
            
        end
        errgauss = [errgauss abs(I-Igauss)];
        errcc    = [errcc abs(I-Icc)];
        
    end
    
	sp = figure(1)
    set(sp,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName',fontname);
    set(sp,'defaultTextFontName', fontname)
    semilogy(nn*M,errcc,['bs' linestyle{k}]);
    grid on
    hold on
    semilogy(nn*M,errgauss,['ro' linestyle{k}])
    grid on
%     if (k==1)
        legendstr{end+1} = ['CC M=' num2str(M)];
        legendstr{end+1} = ['Gauss M=' num2str(M)];
%     else
%         legappend('CC'); 
%         legappend('Gauss');
%     end
    title(testcase)
    % title('Gauss vs. Clenshaw-Curtis quadrature')
end

legend(legendstr,'Location','Best')
xlabel('N_{tot} = M * N')
ylabel('error')
box on

pdfname = testcase;
if (save_pdf==1)
    export_fig(pdfname,'-pdf','-transparent')
end