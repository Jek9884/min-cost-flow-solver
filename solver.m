% TODO passare il vettore di H da ruotare ad apply_old_rotations

% SOLVER Solve a LS problem with GMRES
%
% input:
%   E -> adjacency matrix of a graph without the last row
%   D -> diagonal matrix
%   b -> b vector of the linear system   
%   reorth_flag -> flag to compute the reorthogonalization
%   P -> preconditioning matrix
%
% output:
%   x -> solution
%
function [x] = solver(E, D, b, reorth_flag, P)

    % check if preconditioning is required
    if nargin > 4
        precond_flag = true;
    else
        precond_flag = false;
    end

    % check if reorthogonalization is required
    if isempty(reorth_flag)
        reorth_flag = false;
    end

    % init variables
    k = 1;
    e = [norm(b);0];
    mat_dim = length(D) + length(E);
    givens_list = cell([mat_dim,1]);
    H = [];
    
    % create the first vector of Q from the vector b
    q = b/norm(b);
    Q = q;

    while 1
        if k>mat_dim
            break
        end
        
        % compute the diagonal and sub-diagonal values of H_k and the new
        % column of Q
        [alpha, beta, Q] = compute_new_values(D, E, Q, reorth_flag);
        
        H(k, k) = alpha;
        % check if the Lucky breakdown happens
        if beta < 1e-10
            disp("LUCKY BREAKDOWN")
            break
        end
        H(k + 1, k) = beta;
        H(k, k + 1) = beta;
        

        % apply the old rotation to the new column of H
        H = apply_old_rotations(k, H, givens_list);
        
        % compute the new rotation w.r.t. the new vector of H to zeros the
        % subdiagonal value
        actual_rotation = get_givens_rotation(H(:,k), k, k+1);
        givens_list{k} = actual_rotation;
        H(k:k+1,k) = actual_rotation'*H(k:k+1,k);
        
        % rotate also the ||b||*e1 vector
        e(k:k+1) = actual_rotation'*e(k:k+1);
       
        k = k+1;
        e = [e;0];
    end

    % solve the LS problem
    y = H(1:k-1,1:k-1)\e(1:k-1,1);
    x = Q(:,1:k-1)*y;

end

% GET_GIVENS_ROTATION Compute the Givens rotation matrix w.r.t. the values
% in the position i and k of the vector v. The rotation matrix will contain
% the coefficients [[c, -s];[s, c]]
%
% input:
%   v -> vector to rotate
%   i -> index of the upper value
%   k -> index of the lower value   
%   
% output:
%   rotation_matrix -> matrix used to zeros the lower value
%
function rotation_matrix = get_givens_rotation(v, i, k)

    rotation_matrix = eye(2); 

    % Coefficients to compute the Givens matrix
    a = v(i);
    b = v(k);

    % Compute the Givens matrix
    denominator = sqrt(a^2 + b^2);
    c = a/denominator;
    s = b/denominator;

    % Save the coefficients into the matrix Q
    rotation_matrix(1,1) = c; rotation_matrix(1,2) = -s;
    rotation_matrix(2,1) = s; rotation_matrix(2,2) = c;
     
end

% APPLY_OLD_ROTATION Rotate the k-th vector of H with all the Givens
% rotation matrices of givens_list 
%
% input:
%   k -> index of the vector of H to apply the rotations
%   H -> matrix
%   givens_list -> list of Givens rotation matrices to apply
%   
% output:
%   H -> edited matrix
%
function H = apply_old_rotations(k, H, givens_list)
    % With k=2 it performs just one rotation
    % With k>2 it performs two rotation rotation
    % With k<2 no rotations are applied
    if k ~= 1
        for i = max(1,k-2):(k-1)
            H(i:i+1,k) = givens_list{i}'*H(i:i+1,k);
        end
    end

end

% COMPUTE_NEW_VALUES It computes the diagonal and sub-diagonal values,
% respectively alpha e beta, and apply two times the reorthogonalization 
% (if required)
%
% input:
%   D -> input diagonal matrix (upper-left block of A)
%   E -> input adjcency matrix (lower-left block of A)
%   Q -> list of Givens rotation matrices to apply
%   reorth_flag -> flag to apply reorthogonalization
%   
% output:
%   H -> edited matrix
%
function [alpha, beta, Q] = compute_new_values(D, E, Q, reorth_flag)
    
    % get last row of q
    q = Q(:, end);
    m = length(D);
    % split q in two parts to compute z
    q1 = q(1:m, :);
    q2 = q(m+1:end, :);
    z = [D*q1 + E'*q2; E*q1];

    alpha = q'*z;
    
    if reorth_flag
        % apply reorthogonalization twice
        z = z - Q*(Q'*z);
        z = z - Q*(Q'*z);
    end

    beta = norm(z);
    q = z/beta;
    % append column q to the matrix Q
    Q = [Q,q];
    
end