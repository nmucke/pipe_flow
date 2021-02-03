function type = check_RK_type(A)
% CHECK_RK_TYPE
% returns type of RK method based on Butcher tableau
% type can be: explicit, esdirk, dirk, esimplicit, implicit
s = length(A);
thres = 1e-14;

if (abs(max(max(abs(triu(A)))))<=thres) % upper triangle including diagonal zero
    type = 'explicit';
elseif (abs(max(max(abs(triu(A,1)))))<=thres && max(abs(diag(A)))>=thres && s>1) % upper triangle zero and diagonal nonzero
    if (max(abs(A(1,:)))<=thres)
        type = 'esdirk'; %dirk with explicit first stage
    else
        type = 'dirk';
    end
else
    if (max(abs(A(1,:)))<=thres)
        type = 'esimplicit'; %implicit with explicit first stage 
    else
        type = 'implicit';
    end
end