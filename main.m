seed = 42;
filename = "graphs/net10_32_1.dmx";
[E, D, b] = utility_read_matrix(filename, seed);

disp(size(E));
starting_point  = b;
threshold       = 1e-10;
reorth_flag     = true;
P = NaN;
tic;
[x, e, r_norm] = our_gmres(D, E, P, b, starting_point, threshold, false);
toc;

disp(r_norm);

