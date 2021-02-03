function u = Trefethen1_ex(t,options)
% du/dt = exp(-lambda*t^2)
u = zeros(length(t),1);

lambda  = options.model.constants(1);
u0      = options.model.u0;
t_start = options.time.t_start;

% f       = @(tt) exp(-lambda*tt.^2);
% for i=1:length(t)
%     if (t(i)>t_start)
%         I       = sum(chebfun(f,[t_start,t(i)]));   
%         u(i)    = u0 + I;
%     end
%    
% end

u = (1/2)*sqrt(pi/lambda)*(erf(t*sqrt(lambda)) - ...
                           erf(t_start*sqrt(lambda)));

end