seed = 42;
filename = "graphs/net10_8_1.dmx";
[E, D, b] = utility_read_matrix(filename, seed);

dim = size(D,1)+size(E,1);

starting_point  = b;
max_iterations  = dim;
threshold       = 1e-10;
reorth_flag     = true;
tic;
[x, e, r_norm] = our_gmres(D, E, b, starting_point, max_iterations, threshold, reorth_flag);
toc;

disp(r_norm);

tic;
[x, e, r_norm] = our_gmres(D, E, b, starting_point, max_iterations, threshold, false);
toc;

disp(r_norm);