seed = 42;

[E, D, b] = utility_read_matrix("graphs/net8_8_1.dmx", seed, false);

dim = size(D, 1) + size(E, 1);
starting_point = b;
threshold = 1e-10;
reorth_flag = true;

A = zeros(dim, dim);
A(1:size(D, 1), 1:size(D, 1)) = diag(D);
A(size(D, 1)+1:end, 1:size(E, 2)) = E;
A(1:size(D, 1), size(E, 2)+1:end) = E';

[S, P] = create_schur_complement(D,E);

P_edit = ichol(sparse(P));

[x,r_norm] = our_gmres_slow(A, NaN, b, starting_point, threshold, reorth_flag);
disp(r_norm);

%P_edit_inv_t = P_edit' \ eye(length(P_edit'));
%P_edit_inv = P_edit \ eye(length(P_edit));

%f = @(v) (P_edit_inv_t * A * P_edit_inv)*v;
%b = P_edit_inv_t*b;

%[x, ~, resrel, n_iter] = minres(A, b, 1e-10, dim, P_edit', P_edit, b);
%disp(resrel)
%disp(n_iter)


%{
format long;
seed = 42;
filename = ["graphs/net8_8_1.dmx",
            "graphs/net8_8_2.dmx",
            "graphs/net8_8_3.dmx",
            "graphs/net8_8_4.dmx",
            "graphs/net8_8_5.dmx",
            "graphs/net8_16_1.dmx",
            "graphs/net8_16_2.dmx",
            "graphs/net8_16_3.dmx",
            "graphs/net8_16_4.dmx",
            "graphs/net8_16_5.dmx",
            "graphs/net10_8_1.dmx",
            "graphs/net10_8_2.dmx",
            "graphs/net10_8_3.dmx",
            "graphs/net10_8_4.dmx",
            "graphs/net10_8_5.dmx",
            "graphs/net16_8_1.dmx",
            "graphs/net16_8_2.dmx",
            "graphs/net16_8_3.dmx",
            "graphs/net16_8_4.dmx",
            "graphs/net16_8_5.dmx",];

for i = 1:length(filename)

    [E, D, b] = utility_read_matrix(filename(i), seed, false);
    dim = size(D, 1) + size(E, 1);

    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    fprintf("Filename %s Det: %f\n", filename(i), det(A))

end
%}

%starting_point  = b;
%threshold       = 1e-6;
%reorth_flag     = true;

