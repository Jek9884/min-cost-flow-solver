addpath("./../")
format long;
seed = 42;
filenames = ["../graphs/net8_8_1.dmx", "../graphs/net10_8_3.dmx"];%, "../graphs/net14_64_1.dmx"];
reorth_flags = [false, true];
threshold       = 1e-10;

for i = 1:length(filenames)
    for j = 1:length(reorth_flags)

        [E, D, b] = utility_read_matrix(filenames(i), seed, true);

        starting_point  = b;
        
        tic;
        [x, r_rel, residuals, break_flag] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flags(j));
        toc;
        
        fprintf("Filename: %s Reorth_flag: %d Res rel: %e\n", filenames(i), reorth_flags(j), r_rel)

        disp("---------------------")

        plot_res(residuals, []);

    end
end

function plot_res(residuals, filename)

    figure;
    semilogy(residuals);
    if ~isempty(filename)
        saveas(gcf, filename);
    end

end 