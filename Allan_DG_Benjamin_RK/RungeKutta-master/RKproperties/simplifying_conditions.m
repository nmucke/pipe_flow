function simplifying_conditions(A_RK,b_RK,c_RK)

s = length(c_RK);

eps = 1e-14;

% B
k=1;
while 1
    if ( abs(sum(b_RK.*(c_RK.^(k-1))) - 1/k)>eps )
        disp(['B(' num2str(k-1) ') holds']);
        break
    end
    k = k+1;
end

% C
for i=1:s % c(1) = 0
    k=1;
    while 1

        if ( abs(sum(A_RK(i,:)'.*(c_RK.^(k-1))) - (1/k)*c_RK(i)^k)>eps )
            % condition doesn't hold, so revert to last k for which it
            % holds
            disp(['C(' num2str(k-1) ') holds for i=' num2str(i)]);
            break
        end
        if ( sum(abs(A_RK(i,:))) < eps && abs(c_RK(i))<eps)
            disp(['C(inf) holds for i=' num2str(i)]);
            break
        end
            
        k = k+1;
    end
end

% D
for j=1:s % c(1) = 0
    k=1;
    while 1
        if ( abs(sum(b_RK.*(c_RK.^(k-1)).*A_RK(:,j)) - (b_RK(j)/k)*(1-c_RK(j)^k))>eps )
            disp(['D(' num2str(k-1) ') holds for j=' num2str(j)]);
            break
        end
        if ( sum(abs(A_RK(:,j))) < eps && (b_RK(j)<eps || abs(c_RK(j)-1)<eps))
            disp(['D(inf) holds for j=' num2str(j)]);
            break
        end
                        
%         keyboard
        k = k+1;
        
    end
end

end