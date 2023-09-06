function [S, P, total_time_S] = create_preconditioner(D, E)
    tic;
        dim = size(D, 1) + size(E, 1);
        D = diag(D);
        P = zeros(dim, dim);
        
        tic;
            S = -E * (D\E');
        creation_time_S = toc;

        P(1:size(D, 1), 1:size(D, 1)) = D;
        P(size(D, 1)+1:end, size(E, 2)+1:end) = -S;

        tic;
            S = ichol(sparse(-S));
        factorization_time = toc;
        total_time_S = creation_time_S+factorization_time;
        
        P = ichol(sparse(P));
    x = toc;
    fprintf("Preconditioner created in %f seconds.\n",x);
end