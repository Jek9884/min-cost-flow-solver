function [x, iter] = lanczos_rotation(A, b, M, tol, max_iter)
% Lanczos method with given rotation
% A: coefficient matrix
% b: right-hand side vector
% M: preconditioning matrix
% tol: tolerance for convergence
% max_iter: maximum number of iterations

n = size(A, 1); % Size of the system
x = zeros(n, 1); % Initial guess
r0 = b - A * x; % Initial residual
z0 = M \ r0; % Apply preconditioner to initial residual
v = zeros(n, max_iter+1); % Krylov subspace basis vectors
alpha = zeros(max_iter, 1); % Diagonal elements of the tridiagonal matrix
beta = zeros(max_iter+1, 1); % Subdiagonal elements of the tridiagonal matrix
Q = zeros(n, max_iter+1); % Orthogonal basis vectors
Q(:, 1) = z0 / norm(z0);
iter = 0;

r_norm = norm(r0);
tol = tol * r_norm; % Residual tolerance

while r_norm > tol && iter < max_iter
    iter = iter + 1;
    
    % Arnoldi process
    w = A * Q(:, iter);
    alpha(iter) = Q(:, iter)' * w;
    w = w - alpha(iter) * Q(:, iter);
    beta(iter+1) = norm(w);
    Q(:, iter+1) = w / beta(iter+1);

    % Solve the rotation problem
    [c, s] = given_rotation(alpha(iter), beta(iter+1));
    
    % Update the solution
    x = x + Q(:, 1:iter) * (c * beta(1:iter));
    r = b - A * x;
    z = M \ r;
    
    % Apply the rotation to the residual
    alpha(iter) = c * alpha(iter) + s * beta(iter+1);
    beta(iter+1) = -s * alpha(iter) + c * beta(iter+1);
    
    % Update the norm of the residual
    r_norm = abs(beta(iter+1));
    
    % Print current iteration information
    fprintf('Iteration %d: Residual = %e\n', iter, r_norm);
end

if r_norm <= tol
    fprintf('Lanczos method converged in %d iterations.\n', iter);
else
    fprintf('Lanczos method reached maximum number of iterations without convergence.\n');
end

end

function [c, s] = given_rotation(a, b)
% Given rotation for Lanczos method
% a, b: elements of the tridiagonal matrix

if b == 0
    c = 1;
    s = 0;
elseif abs(b) > abs(a)
    tau = -a / b;
    s = 1 / sqrt(1 + tau^2);
    c = s * tau;
else
    tau = -b / a;
    c = 1 / sqrt(1 + tau^2);
    s = c * tau;
end

end