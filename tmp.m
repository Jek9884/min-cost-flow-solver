path_to_root = "";
experiment_title = "exp_5";
addpath(path_to_root)
format long;
seed = 42;
filenames = ["graphs/net8_8_3.dmx", "graphs/net10_8_3.dmx", "graphs/net12_8_3.dmx"];
threshold = 1e-10;
debug = false;

for i = 1:length(filenames)
    filename = filenames(i);
    [E, ~, b] = utility_read_matrix(path_to_root+filename, seed, debug);
    
    D = ones(size(E,2), 1);
    
    tic;
    [S, P, total_time_S] = create_preconditioner(D,E); 
    P = ichol(sparse(P));
    toc;
end

