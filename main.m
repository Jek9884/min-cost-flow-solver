seed = 42;

[E, D, b] = utility_read_matrix("graphs/net8_8_3.dmx", seed, false);

dim = size(D, 1) + size(E, 1);
starting_point = b;
threshold = 1e-10;
reorth_flag = true;
debug = false;

A = zeros(dim, dim);
A(1:size(D, 1), 1:size(D, 1)) = diag(D);
A(size(D, 1)+1:end, 1:size(E, 2)) = E;
A(1:size(D, 1), size(E, 2)+1:end) = E';

[S, P] = create_preconditioner(D,E);

% ==================== 
trials = 1;

disp("VERSIONE SLOW SENZA PRECONDITIONING")
total_time = 0;
for trial=1:trials
    tic;
    [x,r_rel, k] = our_gmres_slow(A, NaN, b, starting_point, threshold, reorth_flag);
    trial_time = toc;
    total_time = total_time + trial_time;
    fprintf("Trial: %d | Res. Rel: %e | Iter: %d | Trial time : %f\n", trial, r_rel, k, trial_time );
end
fprintf(">>>> Total mean time: %f\n\n", total_time/trials);

disp("VERSIONE SLOW CON PRECONDITIONING")
total_time = 0;
for trial=1:trials
    tic;
    [x,r_rel, k] = our_gmres_slow(A, P, b, starting_point, threshold, reorth_flag);
    trial_time = toc;
    total_time = total_time + trial_time;
    fprintf("Trial: %d | Res. Rel: %e | Iter: %d | Trial time : %f\n", trial, r_rel, k, trial_time );
end
fprintf(">>>> Total mean time: %f\n\n", total_time/trials);

disp("VERSIONE FAST SENZA PRECONDITIONING")
total_time = 0;
for trial=1:trials
    tic;
    [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag, debug);
    trial_time = toc;
    total_time = total_time + trial_time;
    fprintf("Trial: %d | Res. Rel: %e | Iter: %d | Trial time : %f\n", trial, r_rel, k, trial_time );
end
fprintf(">>>> Total mean time: %f\n\n", total_time/trials);


disp("VERSIONE FAST CON PRECONDITIONING")
total_time = 0;
for trial=1:trials
    tic;
    [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, S, b, starting_point, threshold, reorth_flag, false);
    trial_time = toc;
    total_time = total_time + trial_time;
    fprintf("Trial: %d | Res. Rel: %e | Iter: %d | Trial time : %f\n", trial, r_rel, k, trial_time);
end
fprintf(">>>> Total mean time: %f\n\n", total_time/trials);