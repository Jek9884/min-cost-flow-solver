path_to_root = "../../../";
experiment_title = "exp_1";
addpath(path_to_root)
format long;
seed = 42;
filenames       = [path_to_root+"graphs/net8_8_3.dmx" ];%"../graphs/net10_8_3.dmx", "../graphs/net14_64_1.dmx"];
reorth_flags    = [false, true];
threshold       = 1e-10;
debug           = false;

file_path = experiment_title+"_results.csv";
fileID = fopen(file_path, 'w');
fprintf(fileID, "file_name;cond;det;reorth;relative residual;number of iterations;time\n");

for i = 1:length(filenames)
    for j = 1:length(reorth_flags)

        [E, D, b] = utility_read_matrix(filenames(i), seed, debug);

        [c,d] = calculate_det_and_cond(D,E);

        starting_point = b;
        
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flags(j), debug);
        execution_time = toc;

        string_list = split(filenames(i), "/");
        name = string_list(end);
        tmp = split(name, '.');
        name = tmp(1);
        if reorth_flags(j)
            reorth_title = "with_reorth";
        else
            reorth_title = "without_reorth";
        end
        plot_file_name = experiment_title+"_"+name+"_"+reorth_title+".png";
        plot_res(residuals, plot_file_name);


        fprintf(fileID,"%s;%e;%e;%d;%e;%d;%f\n", name, c,d,reorth_flags(j), r_rel, k,execution_time);

    end
end

fclose(fileID);

function plot_res(residuals, filename)

    figure;
    semilogy(residuals);
    if ~isempty(filename)
        saveas(gcf, filename);
    end

end 


function [c,d] = calculate_det_and_cond(D,E)
    dim = size(D, 1) + size(E, 1);

    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    c = cond(A);
    d = det(A);
end
