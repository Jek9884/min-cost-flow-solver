seed = 42;
filename = "graphs/net8_8_3.dmx";
[E, D, b] = utility_read_matrix(filename, seed);

disp(size(E));
starting_point  = b;
threshold       = 1e-10;
reorth_flag     = true;

%tic;
%[x, r_norm] = our_gmres(D, E, NaN, b, starting_point, threshold, false);
%toc;
%disp(r_norm);

S = create_schur_complement(D,E);
%disp(size(S));

%tic;
%[x, r_norm] = our_gmres(D, E, S, b, starting_point, threshold, false);
%toc;
%disp(r_norm);

%S = create_schur_complement(D,E);



