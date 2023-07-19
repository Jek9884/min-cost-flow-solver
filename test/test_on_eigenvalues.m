addpath("./../")
format long;
seed = 42;
filenames = ["../graphs/net8_8_1.dmx"];%, "../graphs/net10_8_3.dmx", "../graphs/net10_32_1.dmx"];
reorth_flag = true;
threshold = 1e-10;

for i = 1:length(filenames)
    
    [E, D, b] = utility_read_matrix(filenames(i), seed, true);

    starting_point  = b;
    
    %Build the full matrix
    dim = size(D, 1) + size(E, 1);
    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';

    spy(eig(A));

    sim_thr = 1e-4;
    n_dist_eig = uniquetol(eig(A), sim_thr);
    fprintf("Number of distinct eigenvalues with a threshold %e: %d\n", sim_thr, length(n_dist_eig));

    tic;
    [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag);
    toc;
    
    fprintf("Filename: %s Reorth_flag: %d Res rel: %e Num iter: %d\n", filenames(i), reorth_flag, r_rel, k)

    disp("---------------------")

end