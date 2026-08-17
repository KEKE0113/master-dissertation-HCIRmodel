function Phat = CIR_Bond_Fixed(r0, lambda, theta, eta, T, N, M)
    dt = T/N;
    r = r0*ones(M,1);
    intR = zeros(M,1);            % ∫ r dt 的累积
    for i = 1:N
        r_pos = max(r,0);
        dW = sqrt(dt)*randn(M,1);
        r_new = r + lambda*(theta - r_pos)*dt + eta*sqrt(r_pos).*dW;
        intR = intR + 0.5*(r + r_new)*dt;   % 梯形积分，比左端点准
        r = r_new;
    end
    Phat = mean(exp(-intR));
end

