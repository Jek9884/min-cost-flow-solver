function [x,H] = original(A, b)
   
    n = 1;
    e = [norm(b);0];
    m = length(A);
    Q_Givens_list = cell([m,1]);
    H = [];

    qq = b/norm(b);
    Q = qq;
    r = A * Q(:, 1);
    
    H(1, 1) = Q(:, 1)' * r; % alpha_1, 1st diagonal element of the tridiagonal matrix
    r = r - H(1, 1) * Q(:, 1); % Remove component in the direction of the first Lanczos vector

    beta = norm(r);
    H(2, 1) = beta; % beta_1, 1st off-diagonal element of the tridiagonal matrix
    H(1, 2) = beta;

    while 1
        if n>m
            break
        end

        z = A*qq;
        alpha = qq'*z;

        %z = z - Q*(Q'*z);
        %z = z - Q*(Q'*z);
            
        beta = norm(z);
        qq = z/beta;
        Q = [Q,qq];

        H(n, n) = alpha;
        if abs(beta) < 1e1
            disp("LUCKY BREAKDOWN")
            break
        end
        H(n + 1, n) = beta;
        H(n, n + 1) = beta;
        

       %QR on H using Givens QR
        if n ~= 1 %Dalla seconda iterazione in poi
            if n>2
               i = n-2;
               H(i:i+1,n) = Q_Givens_list{i}'*H(i:i+1,n);

               i = n-1;
               H(i:i+1,n) = Q_Givens_list{i}'*H(i:i+1,n);
            else
                for i = 1:(n-1)
                    H(i:i+1,n) = Q_Givens_list{i}'*H(i:i+1,n);
                end
            end
        end
        %disp(H);
        Q_Givens = Givens_QR(H(:,n), n, n+1);
        Q_Givens_list{n} = Q_Givens;
        H(n:n+1,n) = Q_Givens'*H(n:n+1,n);
        
        e(n:n+1) = Q_Givens'*e(n:n+1);
       
        n = n+1;
        e = [e;0];
    end

    y = H(1:n-1,1:n-1)\e(1:n-1,1);
    x = Q(:,1:n-1)*y;
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
