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

P = create_schur_complement(D,E);

disp(cond(A));
disp(cond(P\eye(size(P))));
disp(cond(P\A));

%tic;
%[x, r_norm] = our_gmres_slow(A, P, b, starting_point, threshold, false);
%toc;
%disp(r_norm);

%tic;
%[x, r_norm] = our_gmres(D, E, P, b, starting_point, threshold, false);
%toc;
%disp(r_norm);

%tic;
%[x, flag, relres, iter, resvec] = gmres_precond(P, A, b, false, dim, threshold);
%toc;
%disp(norm(b-A*x) / norm(b));

%tic;
%x2 = gmres(A,b, false, threshold, dim, P);
%x2 = lanczos_rotation(A,b,eye(2304,2304),threshold,dim);
%toc;
%r = (b - A*x2);
%disp(norm(r));


%S = create_schur_complement(D,E);



