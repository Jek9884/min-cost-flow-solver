function [x] = solver(E, D, b, varargin)
    
    k = 1;
    e = [norm(b);0];
    m = length(D);
    n = length(E); 
    mat_dim = m + n;
    Q_Givens_list = cell([mat_dim,1]);
    H = [];

    qq = b/norm(b);
    Q = qq;
    qq1 = qq(1:m, :);
    qq2 = qq(m+1:end, :);
    r = [D*qq1 + E'*qq2; E*qq1];
    
    % alpha_1, 1st diagonal element of the tridiagonal matrix
    H(1, 1) = Q(:, 1)' * r;
    % Remove component in the direction of the first Lanczos vector
    r = r - H(1, 1) * Q(:, 1); 

    beta = norm(r);
    % beta_1, 1st off-diagonal element of the tridiagonal matrix
    H(2, 1) = beta; 
    H(1, 2) = beta;

    while 1
        if k>mat_dim
            break
        end
        
        qq1 = qq(1:m, :);
        qq2 = qq(m+1:end, :);
        z = [D*qq1 + E'*qq2; E*qq1];

        alpha = qq'*z;

        %z = z - Q*(Q'*z);
        %z = z - Q*(Q'*z);

        beta = norm(z);
        qq = z/beta;
        Q = [Q,qq];

        H(k + 1, k) = beta;
        H(k, k + 1) = beta;
        H(k, k) = alpha;

        %QR on H using Givens QR
        if k ~= 1 %Dalla seconda iterazione in poi
            if k>2
               i = k-2;
               H(i:i+1,k) = Q_Givens_list{i}'*H(i:i+1,k);

               i = k-1;
               H(i:i+1,k) = Q_Givens_list{i}'*H(i:i+1,k);
            else
                for i = 1:(k-1)
                    H(i:i+1,k) = Q_Givens_list{i}'*H(i:i+1,k);
                end
            end
        end

        Q_Givens = Givens_QR(H(:,k), k, k+1);
        Q_Givens_list{k} = Q_Givens;
        H(k:k+1,k) = Q_Givens'*H(k:k+1,k);
        
        e(k:k+1) = Q_Givens'*e(k:k+1);
       
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
