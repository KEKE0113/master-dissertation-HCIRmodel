function C = HCIR_ChF_Call(S0, K, T, r0, lambda, theta, eta, ...
                                  v0, kappa, vbar, gamma, rho_xv, rho_xr)
% HCIR_ChF_CALL  Semi-analytical European call price under the H1-CIR model
% (Grzelak & Oosterlee 2011, deterministic approximation), using the
% discounted characteristic function (eqs 4.6-4.9 of the thesis, little-trap
% form) and Gil-Pelaez inversion.
%
%   C   : call price
%   P0T : zero-coupon bond price P(0,T) = phi(0)  (useful for parity checks)
%
% Notation follows the thesis / G-O convention:
%   CIR rate:      (lambda, theta, eta),  Heston variance: (kappa, vbar, gamma)
%   correlations:  rho_xv, rho_xr  (rho_vr = 0 imposed)

    k    = log(K);
    P0T  = real(chf(0));           % phi(0) = E[exp(-int r)] = ZCB price
    phmi = real(chf(-1i));         % phi(-i) = S0 (martingale check)

    i1 = @(u) real( exp(-1i*u*k) .* chf(u - 1i) ./ (1i*u*phmi) );
    i2 = @(u) real( exp(-1i*u*k) .* chf(u)      ./ (1i*u*P0T ) );

    I1 = integral(@(u) arrayfun(i1, u), 1e-10, 200, ...
                  'AbsTol', 1e-10, 'RelTol', 1e-8);
    I2 = integral(@(u) arrayfun(i2, u), 1e-10, 200, ...
                  'AbsTol', 1e-10, 'RelTol', 1e-8);

    P1 = 0.5 + I1/pi;              % Q^S  probability  (cf. Phi(d1) in BS)
    P2 = 0.5 + I2/pi;              % T-forward measure (cf. Phi(d2) in BS)

    C  = S0*P1 - K*P0T*P2;

    % ---------- nested: discounted ChF, eqs (4.6)-(4.9) ----------
    function phi = chf(u)
        if eta == 0
            Cf = @(s) (1i*u - 1) .* (1 - exp(-lambda*s)) / lambda;
        else
        D1 = sqrt(lambda^2 + 2*eta^2*(1 - 1i*u));
        G1 = (lambda - D1) / (lambda + D1);
        Cf = @(s) (lambda - D1) .* (1 - exp(-D1*s)) ...
                  ./ (eta^2   * (1 - G1*exp(-D1*s)));
        end
        D2 = sqrt((gamma*rho_xv*1i*u - kappa)^2 - (1i*u - 1)*1i*u*gamma^2);
        G2 = (kappa - gamma*rho_xv*1i*u - D2) / (kappa - gamma*rho_xv*1i*u + D2);

        Df = @(s) (kappa - gamma*rho_xv*1i*u - D2) .* (1 - exp(-D2*s)) ...
                  ./ (gamma^2 * (1 - G2*exp(-D2*s)));
        
        % A(u,T): trapezoidal rule on s in [0, T]
        Ns = 400;
        s  = linspace(0, T, Ns+1);
        Ev = Esqrt(T - s, v0, kappa,  vbar,  gamma);   % E[sqrt(v(T-s))]
        Er = Esqrt(T - s, r0, lambda, theta, eta);     % E[sqrt(r(T-s))]
        f  = kappa*vbar*Df(s) + lambda*theta*Cf(s) ...
             + rho_xr*eta*1i*u .* Ev .* Er .* Cf(s);
        A  = trapz(s, f);

        phi = exp(A + 1i*u*log(S0) + Cf(T)*r0 + Df(T)*v0);
    end
end

% ---------- E[sqrt(x_t)] for a CIR-type process (G-O approximation) ----------
function Ex = Esqrt(t, x0, k_, xbar, g_)
    Ex = zeros(size(t));
    small = (t < 1e-10);
    Ex(small) = sqrt(x0);
    tt  = t(~small);
    if g_ < 1e-12                          % deterministic limit
        Ex(~small) = sqrt(xbar + (x0 - xbar) .* exp(-k_*tt));
        return
    end
    c   = g_^2 * (1 - exp(-k_*tt)) / (4*k_);
    d   = 4*k_*xbar / g_^2;
    lb  = 4*k_*x0*exp(-k_*tt) ./ (g_^2 * (1 - exp(-k_*tt)));
    Ex(~small) = sqrt( c .* (lb - 1 + d + d./(2*(d + lb))) );
end
