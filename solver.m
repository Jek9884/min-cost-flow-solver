function [x] = solver(E, D, b, reorth_flag, P)
    
    if isempty(P)
        precond_flag = false;
    else
        precond_flag = true;
    end

    if isempty(reorth_flag)
        reorth_flag = false;
    end

    k = 1;
    e = [norm(b);0];
    m = length(D);
    n = length(E); 
    mat_dim = m + n;
    givens_list = cell([mat_dim,1]);
    H = [];

    q = b/norm(b);
    Q = q;

    while 1
        if k>mat_dim
            break
        end

        [alpha, beta, q] = compute_new_values(D, E, q, m, reorth_flag, Q);
        Q = [Q,q];
        
        H(k, k) = alpha;
        if beta < 1e-10
            disp("LUCKY BREAKDOWN")
            break
        end
        H(k + 1, k) = beta;
        H(k, k + 1) = beta;
        

        % QR on H using Givens QR
        H = apply_old_rotations(k, H, givens_list);

        actual_rotation = Givens_QR(H(:,k), k, k+1);
        givens_list{k} = actual_rotation;
        H(k:k+1,k) = actual_rotation'*H(k:k+1,k);
        
        e(k:k+1) = actual_rotation'*e(k:k+1);
       
        k = k+1;
        e = [e;0];
    end
    
    y = H(1:k-1,1:k-1)\e(1:k-1,1);
    x = Q(:,1:k-1)*y;

end

function Q = Givens_QR(v, i, k)
    % Initialization of Q, the Givens matrix, as 2x2 identity matrix
        % it will contain the coefficients [[c, -s];[s, c]]
    Q = eye(2); 

    % Coefficients to compute the Givens matrix
    a = v(i);
    b = v(k);

    % Compute the Givens matrix
    c = a/sqrt(a^2 + b^2);
    s = b/sqrt(a^2 + b^2);

    % Save the coefficients into the matrix Q
    Q(1,1) = c;
    Q(2,2) = c;
    Q(2,1) = s;
    Q(1,2) = -s;
end

function H = apply_old_rotations(k, H, givens_list)
    
    if k ~= 1
        for i = max(1,k-2):(k-1)
            H(i:i+1,k) = givens_list{i}'*H(i:i+1,k);
        end
    end

end

function [alpha, beta, q] = compute_new_values(D, E, q, m, reorth_flag, Q)
    
        q1 = q(1:m, :);
        q2 = q(m+1:end, :);
        z = [D*q1 + E'*q2; E*q1];

        alpha = q'*z;

        if reorth_flag
            z = z - Q*(Q'*z);
            z = z - Q*(Q'*z);
        end

        beta = norm(z);
        q = z/beta;
    
end