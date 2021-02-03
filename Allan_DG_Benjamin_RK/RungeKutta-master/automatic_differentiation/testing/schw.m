function y = schw(x)
% 
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% 
n = length(x);
s = sum(-x.*sin(sqrt(abs(x))));
y = 418.9829*n+s;