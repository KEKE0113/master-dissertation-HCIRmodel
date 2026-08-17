function Phat = CIR_Bond_Adap(r0, lambda, theta, eta, tol, hmax, T, M)
    r = r0*ones(M,1);
    t = zeros(M,1);
    intR = zeros(M,1);
    active = true(M,1);
    while any(active)
        idx = find(active);
        h = min(hmax, tol ./ (lambda^2 * max(abs(theta - r(idx)), eps)));
        h = min(h, T - t(idx));
        dW = sqrt(h).*randn(size(h));
        r_old = r(idx);
        r_pos = max(r_old,0);
        r_new = r_old + lambda*(theta - r_pos).*h + eta*sqrt(r_pos).*dW;
        intR(idx) = intR(idx) + 0.5*(r_old + r_new).*h;   % 梯形积分
        r(idx) = r_new;
        t(idx) = t(idx) + h;
        active(idx) = t(idx) < T - 1e-14;
    end
    Phat = mean(exp(-intR));
end

