function [S, P] = create_preconditioner(D, E)
    tic;
        dim = size(D, 1) + size(E, 1);
        D = diag(D);
        P = zeros(dim, dim);
        
        S = -E * (D\E');

        P(1:size(D, 1), 1:size(D, 1)) = D;
        P(size(D, 1)+1:end, size(E, 2)+1:end) = -S;

        P = ichol(sparse(P));
        S = ichol(sparse(-S));
    x = toc;
    fprintf("S created in %f seconds.\n",x);
end