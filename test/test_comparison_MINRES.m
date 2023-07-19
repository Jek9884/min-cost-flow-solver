addpath("./../")
format long;
seed = 42;
filenames = ["../graphs/net8_8_1.dmx"];%, "../graphs/net10_8_3.dmx", "../graphs/net14_64_1.dmx"];
init_mode = ["random"];%, "identity", "all_diff", "ill-conditioned"];
reorth_flags = [false, true];
threshold = 1e-10;

for i = 1:length(filenames)
    
    [E, ~, b] = utility_read_matrix(filenames(i), seed, true);
    starting_point  = b;

    for t = 1:length(init_mode)
              
        D = init_D(size(E,2), init_mode(t));
        
        residuals = {};

        for j = 1:length(reorth_flags)
            
            tic;
            [~, r_rel, res_vec, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flags(j));
            toc;
    
            residuals{j} = res_vec;
            
            fprintf("Filename: %s Reorth_flag: %d Res rel: %e Num iter: %d\n", filenames(i), reorth_flags(j), r_rel, k)
    
        end

        dim = size(D, 1) + size(E, 1);
        A = zeros(dim, dim);
        A(1:size(D, 1), 1:size(D, 1)) = diag(D);
        A(size(D, 1)+1:end, 1:size(E, 2)) = E;
        A(1:size(D, 1), size(E, 2)+1:end) = E';
        
        tic;
        [x, ~, res_rel, n_iter, res_vec] = gmres(A, b, [], threshold);
        toc;
        
        residuals{3} = res_vec/norm(b);
        
        filename = sprintf("%s_%s.png", filenames(i), init_mode(t));
        plot_res(residuals, filename);

        fprintf("Matlab implementation - Filename: %s Res rel: %e Num iter: %d\n", filenames(i), res_rel, n_iter);
        disp("---------------------")
    end

end

function plot_res(residuals, filename)

    tiledlayout(1,3);
    
    nexttile;
    semilogy(residuals{1});
    title("Our GMRES without regularization");
    nexttile;
    semilogy(residuals{2});
    title("Our GMRES with regularization");
    nexttile;
    semilogy(residuals{3});
    title("Matlab MINRES implementation");
    
    if ~isempty(filename)
        saveas(gcf, filename);
    end

end

function [D] = init_D(dim, mode)

    switch mode
        case "random"
            D = rand(dim, 1);
        case "identity"
            D = ones(dim, 1);
        case "ill-conditioned"
            D = logspace(1, -6, dim)';
        case "all_diff"
            D = 1:0.1:(1 + (dim-1)*0.1);
            D = D';
        otherwise
            disp("Init mode not valid, D initialized with uniform distribution")
            D = rand(dim, 1);
    end

end