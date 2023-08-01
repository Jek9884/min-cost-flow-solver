path_to_root = "../../../";
experiment_title = "exp_3";
addpath(path_to_root)
format long;
seed = 42;
filenames = [path_to_root+"graphs/net8_8_1.dmx"];%,path_to_root+"graphs/net8_8_2.dmx",path_to_root+"graphs/net8_8_3.dmx"];%, ,"graphs/net10_8_3.dmx", "graphs/net12_8_3.dmx", "graphs/net14_8_3.dmx"];
threshold = 1e-10;
debug = false;

file_path = experiment_title+"_results.csv";
fileID = fopen(file_path, 'w');
fprintf(fileID, "file_name;cond;det;our relative residual;our number of iterations;our time;GMRES relative residual;GMRES number of iterations;MINRES time;GMRES relative residual;MINRES number of iterations;MINRES time\n");

for i = 1:length(filenames)
    filename = filenames(i);
    [E, D, b] = utility_read_matrix(filename, seed, debug);
    starting_point  = b;
    
    [A,c,d] = calculate_det_and_cond(D,E);
    residuals = {};

    tic;
    [~, our_r_rel, our_res_vec, break_flag, our_k] = our_gmres(D, E, NaN, b, starting_point, threshold, true, debug);
    our_time = toc;

    residuals{1} = our_res_vec(4:end);

    dim = size(D, 1) + size(E, 1);
    
    tic;
    [~, ~, gmres_r_rel, gmres_n_iter, res_vec] = gmres(A, b, [], threshold, dim);
    gmres_time = toc;
    gmres_k = gmres_n_iter(2);

    residuals{2} = res_vec(4:end)/norm(b);

    tic;
    [~, ~, minres_r_rel, minres_n_iter, res_vec] = minres(A, b, threshold, dim);
    minres_time = toc;

    residuals{3} = res_vec(4:end)/norm(b);
    

    string_list = split(filename, "/");
    name = string_list(end);
    tmp = split(name, '.');
    name = tmp(1);
    plot_file_name = experiment_title+"_"+name+"_"+"comparison_with_native_methods.png";
    plot_res(residuals, plot_file_name);
   
    fprintf(fileID,"%s;%e;%e;%e;%d;%f;%e;%d;%f;%e;%d;%f;\n", name,c,d, our_r_rel, our_k, our_time, ...
                    gmres_r_rel, gmres_k, gmres_time, ...
                    minres_r_rel, minres_n_iter, minres_time);
end

fclose(fileID);

function plot_res(residuals, filename)
    colors = ["#0072BD","#D95319","#EDB120"];
    figure;
    semilogy(cell2mat(residuals(1)), 'LineWidth',2);
    hold on;
    for i =2:numel(residuals)
       semilogy(cell2mat(residuals(i)), 'LineWidth',2);
    end
    legend(["Our method","MATLAB GMRES","MATLAB MINRES"]);
    hold off;
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end

function [A, c,d] = calculate_det_and_cond(D,E)
    dim = size(D, 1) + size(E, 1);

    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    c = cond(A);
    d = det(A);
end

