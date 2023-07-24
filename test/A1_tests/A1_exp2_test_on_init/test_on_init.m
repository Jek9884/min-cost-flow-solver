path_to_root = "../../../";
experiment_title = "exp_2";
addpath(path_to_root)
format long;
seed = 42;
filenames = [path_to_root+"graphs/net8_8_3.dmx"];
reorth_flag = true;
init_mode = ["random", "identity", "ill-conditioned", "all_diff"];
threshold = 1e-10;
debug = false;
for j = 1:length(filenames)
    filename = filenames(j);
    for i = 1:length(init_mode)
        
        [E, ~, b] = utility_read_matrix(filename, seed, debug);
    
        D = init_D(size(E,2), init_mode(i));
    
        starting_point  = b;
        
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag, debug);
        toc;
        
        fprintf("Filename: %s\t| Init mode: %s\t| Reorth_flag: %d\t| Res rel: %e\t| Num iter: %d\n", filename, init_mode(i), reorth_flag, r_rel, k)
    
        disp("---------------------")
    
        string_list = split(filename, "/");
        name = string_list(end);
        tmp = split(name, '.');
        name = tmp(1);
        plot_file_name = experiment_title+"_"+name+"_"+init_mode(i)+".png";
        plot_res(residuals, plot_file_name);
    end
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