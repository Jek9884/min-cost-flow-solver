addpath("./../")
format long;
seed = 42;
filename = "../graphs/net8_8_3.dmx";
reorth_flag = true;
init_mode = ["random", "identity", "ill-conditioned", "all_diff"];
threshold = 1e-10;

for i = 1:length(init_mode)
    
    [E, ~, b] = utility_read_matrix(filename, seed, true);

    D = init_D(size(E,2), init_mode(i));

    starting_point  = b;
    
    tic;
    [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag);
    toc;
    
    fprintf("Filename: %s Init mode: %s\n Reorth_flag: %d Res rel: %e Num iter: %d\n", filename, init_mode(i), reorth_flag, r_rel, k)

    disp("---------------------")

    %plot_res(residuals, []);

end

function plot_res(residuals, filename)

    figure;
    semilogy(residuals);
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