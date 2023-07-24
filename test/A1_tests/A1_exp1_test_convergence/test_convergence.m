path_to_root = "../../../";
experiment_title = "exp_1";
addpath(path_to_root)
format long;
seed = 42;
filenames = [path_to_root+"graphs/net8_8_3.dmx" ];%"../graphs/net10_8_3.dmx", "../graphs/net14_64_1.dmx"];
reorth_flags    = [false, true];
threshold       = 1e-10;
debug           = false;

for i = 1:length(filenames)
    for j = 1:length(reorth_flags)

        [E, D, b] = utility_read_matrix(filenames(i), seed, debug);

        starting_point = b;
        
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flags(j), debug);
        toc;
        
        fprintf("Filename: %s\t| Reorth_flag: %d\t| Res rel: %e\t| Num iter: %d\n", filenames(i), reorth_flags(j), r_rel, k)

        disp("---------------------")

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

    end
end

function plot_res(residuals, filename)

    figure;
    semilogy(residuals);
    if ~isempty(filename)
        saveas(gcf, filename);
    end

end 