seed = 42;
filename = "graphs/net8_8_1.dmx";
[A, E, D, b] = utility_read_matrix(filename, seed);

disp(length(b))
disp(length(A))

% Compute the largest k eigenvalues/eigenvectors
[V, D] = eigs(A, 5, 'largestabs');

% Construct a matrix A with lucky breakdown
A = V * D * V';

if isequal(A, A')
    disp("SI LO E'")
end
%{
disp("MATLAB MINRES")
tic;
[x] = minres(A,b);
toc;
res = b - A*x;
disp(norm(res))

disp("Custom GMRES with D and E")
tic;
[x] = solver(E, D, b, false);
toc;
res = b - A*x;
disp(norm(res))
%}
disp("Custom GMRES with A")
tic;
[x_og] = original(A,b);
toc;
res_og = b - A*x_og;
disp(norm(res_og))
