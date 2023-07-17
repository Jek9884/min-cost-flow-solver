seed = 42;
filename = "graphs/net8_8_1.dmx";
[E, D, b] = utility_read_matrix(filename, seed);

disp(size(E));
starting_point  = b;
threshold       = 1e-10;
reorth_flag     = true;

dim = size(D, 1) + size(E, 1);

A = zeros(dim, dim);

A(1:size(D, 1), 1:size(D, 1)) = diag(D);
A(size(D, 1)+1:end, 1:size(E, 2)) = E;
A(1:size(D, 1), size(E, 2)+1:end) = E';

[x,r_norm] = our_gmres_slow(A, b, starting_point, threshold, reorth_flag);
disp(r_norm);

