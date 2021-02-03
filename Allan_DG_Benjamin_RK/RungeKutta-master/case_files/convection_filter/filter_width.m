function [ H, Hx ] = filter_width( x, H0, type )
%FILTER_WIDTH Get filter size and its derivative
%   symmetric filtering of a signal from x-H(x) ... x+H(x), with effective
%   filter size 2*H and derivative Hx = dH(x)/dx

switch type
    
    case 'cos'
        %% nonuniform filter width, such that we have 'refinement' at the boundaries
        H     = (-0.5*cos(x)+1)*H0;
        Hx    = (0.5*sin(x))*H0;
        
    case 'uni'
        %% uniform:
        H  = H0*ones(length(x),1);
        Hx = zeros(length(x),1);
        
    otherwise
        error('filter type unknown');
end

end

