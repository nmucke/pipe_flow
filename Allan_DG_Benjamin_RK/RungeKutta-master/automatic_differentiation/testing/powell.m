function y = powell(x)  
% 
% Powell function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 24.
% 
n = length(x);
m = n;
for i = 1:m/4
   fvec(4*i-3) = x(4*i-3)+10*(x(4*i-2));
   fvec(4*i-2) = sqrt(5)*(x(4*i-1)-x(4*i));
   fvec(4*i-1) = (x(4*i-2)-2*(x(4*i-1)))^2;
   fvec(4*i)   = sqrt(10)*(x(4*i-3)-x(4*i))^2;
end;
fvec = fvec';
y = norm(fvec)^2;