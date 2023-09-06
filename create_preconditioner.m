function [S, P, creation_time_S] = create_preconditioner(D, E)
    tic;
        dim = size(D, 1) + size(E, 1);
        D = diag(D);
        P = zeros(dim, dim);
        
        tic;
            S = -E * (D\E');
        creation_time_S = toc;

        P(1:size(D, 1), 1:size(D, 1)) = D;
        P(size(D, 1)+1:end, size(E, 2)+1:end) = -S;
        
        P = ichol(sparse(P));
    total_precond_time = toc;
    fprintf("S matrix created in %f seconds.\n",creation_time_S);
end