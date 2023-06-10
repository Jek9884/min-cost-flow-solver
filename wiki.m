seed = 42;
filename = "graphs/net8_8_1.dmx";
[A, E, D, b] = utility_read_matrix(filename, seed);
tic;
[x, e] = gmres( A, b, b, size(A,1) , 1e-10);
toc;
disp(norm(A*x-b)/norm(b));

function [x,e] = gmres( A, b, starting_point, max_iterations, threshold)
  m = max_iterations;

  r = (b - (A * starting_point));

  b_norm = norm(b);

  sn = zeros(m, 1);
  cs = zeros(m, 1);
  
  r_norm = norm(r);
  Q(:,1) = r / r_norm;

  e1 = zeros(m+1, 1);
  e1(1) = 1;
  e = r_norm * e1;
  
  for k = 1:m

    % run arnoldi
    [H(1:k+1, k), Q(:, k+1)] = arnoldi(A, Q, k);
    
    % eliminate the last element in H ith row and update the rotation matrix
    [H(1:k+1, k), cs(k), sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
   
    temp   =  cs(k) * e(k) + sn(k) * e(k + 1);
    e(k + 1) = -sn(k) * e(k) + cs(k) * e(k + 1);
    e(k)   = temp;
    
    y = H(1:k, 1:k) \ e(1:k);
    x = starting_point + Q(:, 1:k) * y;
    r = norm(b - (A * x))/b_norm;

    if abs(r)<threshold
        fprintf("Terminated in %d iterations\n", k);
        break;
    end
  end

  y = H(1:k, 1:k) \ e(1:k);
  x = starting_point + Q(:, 1:k) * y;
end

function [h, z] = arnoldi(A, Q, k)
  z = A*Q(:,k);  
  for i = max(1,k-2):k    
    h(i) = z' * Q(:, i);
    z = z - h(i) * Q(:, i);
  end
  h(k + 1) = norm(z);
  z = z / h(k + 1);
end

function [h, cs_k, sn_k] = apply_givens_rotation(h, cs, sn, k)
  % apply for ith column
  for i = max(1,k-2):k-1
    temp   =  cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
  end

  % update the next sin cos values for rotation
  [cs_k, sn_k] = givens_rotation(h(k), h(k + 1));

  % eliminate H(i + 1, i)
  h(k) = cs_k * h(k) + sn_k * h(k + 1);
  h(k + 1) = 0.0;
end

%%----Calculate the Givens rotation matrix----%%
function [cs, sn] = givens_rotation(v1, v2)
%  if (v1 == 0)
%    cs = 0;
%    sn = 1;
%  else
    t = sqrt(v1^2 + v2^2);
%    cs = abs(v1) / t;
%    sn = cs * v2 / v1;
    cs = v1 / t;  % see http://www.netlib.org/eispack/comqr.f
    sn = v2 / t;
%  end
end