seed = 42;
filename = "graphs/net8_8_1.dmx";
[E, D, b] = utility_read_matrix(filename, seed);

disp(size(E));
starting_point  = b;
threshold       = 1e-10;
reorth_flag     = true;
P = NaN;

tic;
[x, r_norm] = our_gmres(D, E, P, b, starting_point, threshold, false);
toc;
disp(r_norm);


tic;
[x, r_norm] = our_gmres(D, E, P, b, starting_point, threshold, true);
toc;
disp(r_norm);


P = create_preconditioner(D,E);



