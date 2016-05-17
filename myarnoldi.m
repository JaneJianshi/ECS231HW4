function [V, H, iter] = myarnoldi(A, v0, m, reorth)

n = size(A, 2);
 v_1 = v0 / norm(v0);
 V(1:n,1:(m+1)) = zeros(n,m+1);
    V(:, 1) = v_1;
    H(1:(m+1),1:m) = zeros(m+1,m);
    for j = 1: m
        w = A * V(:, j);
        for i = 1:j
            H(i, j) = V(:, i)' * w; 
            w = w - H(i,j) * V(:, i); 
        end
        
        if reorth == 1
            for i = 1:j
                temp =  V(:, i)' * w; 
                w = w - temp * V(:, i);
                H(i, j) = H(i, j) + temp; 
            end
        end
        

        H(j+1, j) = norm(w); 
        V(:, j+1) = w/H(j+1, j);
        iter = j;
        
        if (H(j+1, j) == 0)
            return; 
        end
    end
    