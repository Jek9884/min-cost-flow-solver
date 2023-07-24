path_to_root = "../../../";
experiment_title = "exp_4";
addpath(path_to_root)
format long;
seed = 42;
filenames = [path_to_root+"graphs/net8_8_3.dmx"];%, "../graphs/net10_8_3.dmx", "../graphs/net14_64_1.dmx"];
threshold = 1e-10;
debug = false;

for i = 1:length(filenames)
    disp("==========================");
    disp(filenames(i))
    disp("==========================");

    [E, D, b] = utility_read_matrix(filenames(i), seed, debug);
    starting_point  = b;
        
    residuals = {};

    tic;
    [~, r_rel, res_vec, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, true, debug);
    time = toc;

    residuals{1} = res_vec;
    fprintf("[OUR IMPLEMENTATION] \tRes rel: %e\t| Num iter: %d\t| Time: %f\n", r_rel, k, time)


    dim = size(D, 1) + size(E, 1);
    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    tic;
    [x, ~, res_rel, n_iter, res_vec] = gmres(A, b, [], threshold, dim);
    time = toc;

    fprintf("[MATLAB GMRES] \t\tRes rel: %e\t| Num iter: %d\t| Time: %f\n", res_rel, n_iter(2),time);
    residuals{2} = res_vec/norm(b);

    tic;
    [x, ~, res_rel, n_iter, res_vec] = minres(A, b, threshold, dim);
    time = toc;

    fprintf("[MATLAB MINRES] \tRes rel: %e\t| Num iter: %d\t| Time: %f\n", res_rel, n_iter, time);
    residuals{3} = res_vec/norm(b);
    

    string_list = split(filename, "/");
    name = string_list(end);
    tmp = split(name, '.');
    name = tmp(1);
    plot_file_name = experiment_title+"_"+name+"_"+"comparison_with_native_methods.png";
    plot_res(residuals, filename);
   
end

function plot_res(residuals, filename)
    tiledlayout(1,3);
 
    nexttile;
    semilogy(residuals{1});
    title("Our implementation");
    nexttile;
    semilogy(residuals{2});
    title("Native GMRES");
    nexttile;
    semilogy(residuals{3});
    title("Native MINRES");
    
    if ~isempty(filename)
        saveas(gcf, filename);
    end

end
