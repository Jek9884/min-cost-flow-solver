function [x, flag, relres, iter, resvec] = gmres_precond(P, A, b, restart, maxit, tol)

    n = size(A, 1);
    m = n;% min(restart, n);

    x = zeros(n, 1);
    r0 = b - A * x;
    z0 = P \ r0;

    beta = norm(z0);
    v = z0 / beta;

    H = zeros(m+1, m);
    cs = zeros(m, 1);
    sn = zeros(m, 1);
    e1 = zeros(m+1, 1);
    e1(1) = 1;

    resvec = zeros(maxit, 1);
    resvec(1) = norm(r0);
    relres = resvec(1) / norm(b);

    iter = 1;
    flag = 0;

    while (iter <= maxit) && (relres > tol)
        v = v / norm(v);
        w = P \ (A * (P' * v));

        for j = 1:m
            fprintf("it: %d\n", j);
            H(j, j) = w' * v;
            w = w - H(j, j) * v;
            for i = 1:j-1
                H(i, j) = w' * v;
                w = w - H(i, j) * v;
            end
            H(j+1, j) = norm(w);
            if H(j+1, j) == 0
                disp("break 1");
                break;
            end
            v = w / H(j+1, j);
            for i = 1:j-1
                temp = cs(i) * H(i, j) + sn(i) * H(i+1, j);
                H(i+1, j) = -sn(i) * H(i, j) + cs(i) * H(i+1, j);
                H(i, j) = temp;
            end
            [cs(j), sn(j)] = givens(H(j, j), H(j+1, j));
            H(j, j) = cs(j) * H(j, j) + sn(j) * H(j+1, j);
            H(j+1, j) = 0;
            e1(j+1) = -sn(j) * e1(j);
            e1(j) = cs(j) * e1(j);
            resvec(iter+1) = norm(A*x-b);
            relres = resvec(iter+1) / norm(b);
            disp(relres);
            iter = iter + 1;
            if (relres <= tol)
                disp("break 2");
                break;
            end
        end
        if (relres <= tol)
            break;
        end
        y = H(1:m, 1:m) \ e1(1:m);
        x = x + v(1:m) * y;
        r = b - A * x;
        disp(r);
        z = P \ r;
        beta = norm(z);
        v = z / beta;
    end

    if (relres > tol)
        flag = 1;
    end

end

function [c, s] = givens(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            temp = a / b;
            s = 1 / sqrt(1 + temp^2);
            c = temp * s;
        else
            temp = b / a;
            c = 1 / sqrt(1 + temp^2);
            s = temp * c;
        end
    end
end
