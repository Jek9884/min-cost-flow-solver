function [x,e,r_norm] = our_gmres(D, E, b, starting_point, max_iterations, threshold)
    
  patient_tol = 1e-12; 

  m = max_iterations;

  % Here we are calculating r = (b - (A * starting_point)) exploiting 
  % the matrix structure.
  sp1 = starting_point(1:size(D,1));
  sp2 = starting_point(size(D,1)+1:end);
  r = (b - [(D.*sp1)+(E'*sp2); E*sp1]);

  b_norm = norm(b);

  sn = zeros(m, 1);
  cs = zeros(m, 1);
  
  r_norm = norm(r);
  Q(:,1) = r / r_norm;

  e1 = zeros(m+1, 1);
  e1(1) = 1;
  e = r_norm * e1;

  last_residual = -1;
  
  for k = 1:m
    [H(1:k+1, k), Q(:, k+1)] = lanczos(D, E,  Q, k, true);
    
    [H(:,k)] = apply_old_rotations(H(:,k), k, cs, sn);

    [H(1:k+1, k),e, cs(k), sn(k)] = apply_current_rotation(H(1:k+1,k), e, k);
    
    y = H(1:k, 1:k) \ e(1:k);
    x = starting_point + Q(:, 1:k) * y;

    x1 = x(1:size(D,1));
    x2 = x(size(D,1)+1:end);
    r = (b - [(D.*x1)+(E'*x2); E*x1]);
    r_norm = norm(r)/b_norm;

    if abs(r_norm)<threshold
        fprintf("Terminated in %d iterations\n", k);
        break;
    end
    
    if abs(last_residual-r_norm) < patient_tol
        patient = patient+1;
    else
        patient = 0;
    end
    if patient >= 3
        fprintf("Terminated in %d iterations (due to the patient) \n", k);
        break;
    end
    
    last_residual = r_norm;
  end

  y = H(1:k, 1:k) \ e(1:k);
  x = starting_point + Q(:, 1:k) * y;
end

function [h, z] = lanczos(D, E, Q, k, reorth_flag)
  q1 = Q(1:size(D,1), k);
  q2 = Q(size(D,1)+1:end, k);
  z = [(D.*q1)+(E'*q2); E*q1];

  for i = max(1,k-2):k    
    h(i) = z' * Q(:, i);
    z = z - h(i) * Q(:, i);
  end

  if reorth_flag  % apply reorthogonalization twice, if necessary
    z = z - Q*(Q'*z);
    z = z - Q*(Q'*z);
  end

  h(k + 1) = norm(z);
  z = z / h(k + 1);
end

%
% This function applies the previous rotations to the newly computed column 
%
% Input: h - TODO
%        k - TODO
%        cs, sn - coefficients of the previous rotations
%
% Output: h - updated vector 
%
function [h] = apply_old_rotations(h, k, cs, sn)
  for i = max(1,k-2):k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
  end
end

% TODO: commenti
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