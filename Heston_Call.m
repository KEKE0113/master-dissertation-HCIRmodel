function C = Heston_Call(S0,K,r,T,v0,kappa,theta,xi,rho)
%HESTON_CALL Semi-analytical European call price under the Heston model
%   C = S0*P1 - K*exp(-r*T)*P2

x0 = log(S0);
k  = log(K);

% Risk-neutral characteristic function of ln(S_T), little trap form
    function phi = cf(u)
        u  = u(:).';                       % row vector
        iu = 1i*u;
        d  = sqrt((kappa - rho*xi*iu).^2 + xi^2*(iu + u.^2));
        g  = (kappa - rho*xi*iu - d) ./ (kappa - rho*xi*iu + d);
        e  = exp(-d*T);
        A  = iu*(x0 + r*T) ...
           + kappa*theta/xi^2 * ((kappa - rho*xi*iu - d)*T ...
             - 2*log((1 - g.*e)./(1 - g)));
        D  = (kappa - rho*xi*iu - d)/xi^2 .* (1 - e)./(1 - g.*e);
        phi = exp(A + D*v0);
    end

% Gil-Pelaez probabilities
% P2: risk-neutral probability Q(S_T > K)
integrand2 = @(u) real(exp(-1i*u*k).*cf(u)./(1i*u));
P2 = 1/2 + 1/pi * integral(integrand2, 1e-8, 200);
 
% P1: probability under the stock numeraire measure Q^S
% cf under Q^S: cf(u - i)/cf(-i), and cf(-i) = S0*exp(r*T)
integrand1 = @(u) real(exp(-1i*u*k).*cf(u - 1i)./(1i*u*S0*exp(r*T)));
P1 = 1/2 + 1/pi * integral(integrand1, 1e-8, 200);
 
C = S0*P1 - K*exp(-r*T)*P2;

end

