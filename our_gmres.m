%
% This is our implementation of the GMRES algorithm, exploiting the matrix structure. 
%

%break_flag is 0 if the method reaches the convergence with the threshold, 1 if converges with the patient and 2 otherwise
function [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, P, b, starting_point, threshold, reorth_flag)

  % Checks on the input parameters
  assert(size(D,2)==1, "D must be a vector");
  assert(size(starting_point,2)==1, "starting_point must be a vector");
  assert(islogical(reorth_flag), "reorth_flag must be a boolean value");
  
  residuals = [];

  dim = size(D,1)+size(E,1);

  % Initialization
  patient_tol = 1e-19; 
  patient = 0;
  m = dim;

  b_norm = norm(b);
 
  r = calculate_the_residual_optimized(D, E, P, b, starting_point);

  % Initialization of the vectors
  sn = zeros(m, 1); cs = zeros(m, 1); % Instead of saving the whole rotation matrix, we save only the coefficients.
  e1 = zeros(m, 1);     
  
  r_norm = norm(r);

  Q(:,1) = r / r_norm;

  e1(1) = 1;
  e = r_norm * e1;

  last_residual = -1;
  
  for k = 1:m
    % Step 1: Lanczos algorithm
    [H(1:k+1, k), Q(:, k+1)] = lanczos(D, E, P, Q, k, reorth_flag);

    % Step 2: apply the previous rotations to the newly computed column
    H(:,k) = apply_old_rotations(H(:,k), k, cs, sn);

    % Step 3: apply the current rotation to the newly computed column
    [H(1:k+1, k), e, cs(k), sn(k)] = apply_current_rotation(H(1:k+1,k), e, k);
    
    % Step 4: check the residual
    y = H(1:k, 1:k) \ e(1:k); 
    x = starting_point + Q(:, 1:k) * y;

    r = calculate_the_residual_optimized(D, E, P, b, x);
    r_rel = norm(r)/b_norm;

    residuals(end+1) = r_rel;

    if abs(r_rel) < threshold
        break_flag = 0;
        break;
    end
    
    % ======== PATIENT ========
    if abs(last_residual - r_rel) < patient_tol
        patient = patient+1;
    else
        patient = 0;
    end
    if patient >= 3
        break_flag = 1;
        break;
    end
    
    last_residual = r_rel;

    %fprintf("Iter: %d Res rel: %e\n", k, last_residual);

  end
    
  if k == m
    break_flag = 2;
  end

  y = H(1:k, 1:k) \ e(1:k);             %O( ([k=]m+n)^2 )
  x = starting_point + Q(:, 1:k) * y;   %O( (m+n)^2 )
end

%
% This function computes the k-th column of the Hessenberg matrix H and the k-th column of the matrix Q
%
% Input: D - the original diagonal vector
%        E - the original E matrix
%        S - TODO
%        Q - TODO
%        k - the index of the column to be computed
%        reorth_flag - this flag is used to decide if the reorthogonalization is needed
%
% Output: h - TODO
%         v - TODO
%
function [h, v] = lanczos(D, E, P, Q, k, reorth_flag)
  q1 = Q(1:size(D,1), k);
  q2 = Q(size(D,1)+1:end, k);

  if ~isnan(P)
     P_edit_inv = inv(P);
     m = size(E,2);
     B1  = (P_edit_inv(1:m,1:m)' .* D)  *  P_edit_inv(1:m,1:m);
     B2  = (P_edit_inv(1:m,1:m)'  * E') *  P_edit_inv((m+1):end,(m+1):end);
     B3  = (P_edit_inv((m+1):end,(m+1):end)'*  E)  *  P_edit_inv(1:m,1:m);
     v = [(B1*q1)+(B2*q2); B3*q1];
  else
     v = [(D.*q1)+(E'*q2); E*q1];
  end

  for i = max(1,k-2):k    
    h(i) = v' * Q(:, i);
    v = v - h(i) * Q(:, i);
  end

  if reorth_flag  % apply reorthogonalization twice, if necessary
    v = v - Q*(Q'*v);
    v = v - Q*(Q'*v);
  end

  h(k + 1) = norm(v);
  v = v / h(k + 1);
end

%
% This function applies the previous rotations to the newly computed column 
%
% Input: h - TODO
%        k - the index of the column to be computed
%        cs, sn - coefficients of the previous rotations
%
% Output: h - updated vector 
%
function h = apply_old_rotations(h, k, cs, sn)
  % Rotations are computed only on the last items of the vector, since TODO 
  for i = max(1,k-2):k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
  end
end

% This function applies the current rotation to the newly computed column
%
% Input: h - the vector to be rotated
%        e - the canonical vector to be rotated
%        k - the index of the column to be computed
% 
% Output: h - updated vector
%         e - updated vector
%         cs_k, sn_k - coefficients of the current rotation
% 
function [h, e, cs_k, sn_k] = apply_current_rotation(h,e,k)
  [cs_k, sn_k] = get_givens_rotation(h,k,k+1);

  tmp       = cs_k * h(k)  + sn_k * h(k + 1);
  h(k + 1)  = -sn_k * h(k) + cs_k * h(k + 1);
  h(k) = tmp;

  temp     =  cs_k * e(k) + sn_k * e(k + 1);
  e(k + 1) = -sn_k * e(k) + cs_k * e(k + 1);
  e(k)     = temp;
end

% 
% This function computes the coefficients of the Givens rotation matrix
% that zeros the j-th component of the vector v.
%
% Input: v    - vector of size 2
%        i,j  - indices of the vector
% Output: c,s - coefficients of the Givens matrix
% 
function [c,s] = get_givens_rotation(v, i, j) % Coefficients to compute the Givens matrix
    a = v(i);
    b = v(j);

    denominator = sqrt(a^2 + b^2);
    c = a/denominator;
    s = b/denominator;
end

%
% By using the matrix structure, the calculation of the residual can be optimized.
%
% Input: D - the original diagonal vector
%        E - the original E matrix
%        S - TODO
%        b - the original b vector
%        input_vector - the vector on which the residual is computed
%
% Output: r - the residual
%
function r = calculate_the_residual_optimized(D, E, P, b, input_vector)
  part_1 = input_vector(1:size(D,1));
  part_2 = input_vector(size(D,1)+1:end);

  if ~isnan(P) 
     P_edit_inv = inv(P);
     m = size(E,2);
     B1  = (P_edit_inv(1:m,1:m)' .* D)  *  P_edit_inv(1:m,1:m);
     B2  = (P_edit_inv(1:m,1:m)'  * E') *  P_edit_inv((m+1):end,(m+1):end);
     B3  = (P_edit_inv((m+1):end,(m+1):end)'*  E)  *  P_edit_inv(1:m,1:m);
    
     b_part_1 = b(1:size(D,1));
     b_part_2 = b(size(D,1)+1:end);

     b = [P_edit_inv(1:m,1:m)*b_part_1; P_edit_inv(m+1:end,m+1:end)*b_part_2];
  
     r = (b - [(B1*part_1)+(B2*part_2); B3*part_1]);
  else
     r = (b - [(D.*part_1)+(E'*part_2); E*part_1]);
  end
end

