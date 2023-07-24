path_to_root = "../../../";
experiment_title = "exp_3";
addpath(path_to_root)
format long;
seed = 42;
filenames = [path_to_root+"graphs/net8_8_3.dmx"];%, "../graphs/net10_8_3.dmx", "../graphs/net10_32_1.dmx"];
reorth_flag = true;
threshold = 1e-10;
debug = false;
distinct_values_on_D = [1,2,5,10,20];

for i = 1:length(filenames)

    disp("==========================");
    disp(filenames(i))
    disp("==========================");

    for j = 1:length(distinct_values_on_D)
        actual_distinct_values_on_D = distinct_values_on_D(j);
        [E, ~, b] = utility_read_matrix(filenames(i), seed, debug);
        D = generate_D(actual_distinct_values_on_D, size(E, 2));
        
        starting_point  = b;
        
        %Build the full matrix
        dim = size(D, 1) + size(E, 1);
        A = zeros(dim, dim);
        A(1:size(D, 1), 1:size(D, 1)) = diag(D);
        A(size(D, 1)+1:end, 1:size(E, 2)) = E;
        A(1:size(D, 1), size(E, 2)+1:end) = E';
    
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag, debug);
        toc;
        
        fprintf("Res rel: %e | Num iter: %d\t| Distict Values on D:%d \n", r_rel, k, actual_distinct_values_on_D);
    
        disp("---------------------")
    end
end


function [D] = generate_D(distinct_values, dim)
    D = ones(dim, 1);
    if distinct_values == 1
        return;
    elseif distinct_values == 2
        %In this case the first half are all ones,
        D(round((dim/2)+1):end) = 2; % while the second half are all twos.
    else
        for i = 2:distinct_values
            start_index = round(((dim/distinct_values)*(i-1))+1);
            end_index = round((dim/distinct_values)*(i));
           
            D(start_index:end_index) = i;
        end
    end
end