% D[HankelH2[n, x], x]
% 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])
% D[BesselJ[n, x], x]
% 1/2 (HankelH2[-1 + n, x] - HankelH2[1 + n, x])

if strcmp(mode, 'TE Hz')
    ka = r*k0;
    Hz_inc = H_inc_0*exp(-1j*k0*r*cos(theta));
    
    %-% Equ 6.4.19
bn = @(n, ka) -(1j)^(-n)*( ...
  ( 0.5*(besselj(n-1,ka)   - besselj(n+1,ka))   )/ ...
  ( 0.5*(besselh(n-1,2,ka) - besselh(n+1,2,ka)) )  ...
  );

    maxN = 40;
    b_n = zeros(2*maxN+1,1);
    for n=-maxN:maxN
        b_n(n+maxN+1) = bn(n, ka);
    end
    figure(21);
    plot(-maxN:maxN, real(b_n), -maxN:maxN, imag(b_n));
    xlabel('index n');title('b_n');legend({'Re(b_n)','Im(b_n)'})

    b_n_H2 = zeros(2*maxN+1, 1);
    for n=-maxN:maxN
        b_n_H2(n+maxN+1) = b_n(n+maxN+1)*besselh(n,2,ka);
    end
    exp_theta = exp(1j* (-maxN:maxN).*theta);

    % n x Hz z = - Hz t = - Hz phi
    J_ideal = -(H_inc_0 * exp_theta * b_n_H2 + Hz_inc);
end

if strcmp(mode, 'TM Ez')
% N = 500;
% theta = linspace(0, 2*pi, N+1).'; theta(end) = []; % column
    ka = r*k0;
    curlE_inc_phi = cos(theta).*(1j*k0*E_inc_0*exp(-1j*k0*r*cos(theta)));
    
    %-% Equ 6.4.12
an_H2p = @(n, ka) -(1j)^(-n)*( ...
( besselj(n,ka) / besselh(n,2,ka) ) * ...
( 0.5*(besselh(n-1,2,ka) - besselh(n+1,2,ka)) )  ...
);
    maxN = 40;
    a_n_H2p = zeros(2*maxN+1,1);
    for n=-maxN:maxN
        a_n_H2p(n+maxN+1) = an_H2p(n, ka);
    end
    figure(21);
    plot(-maxN:maxN, real(a_n_H2p), -maxN:maxN, imag(a_n_H2p));
    xlabel('index n');title('a_n \cdot H^{(2)}_n(k\rho) | \rho = radius'); legend({'Re()','Im()'})
    
    exp_theta = exp(1j* (-maxN:maxN).*theta);
    curlE_scatt_phi = -exp_theta*(a_n_H2p*k0); 
    % - sign due to 
    % \nabla \times E^{scatt} = \frac{1}{\rho} \diff{E_z}{\phi} \hat{\rho}
    % + -\diff{E_z}{\rho} \hat{\phi}
    
    % (-1j*omega*mu) H = \curl E
    % J = n (rho) x (H_rho, H_phi, 0) = \rho x H_phi \phi = H_phi z
    J_ideal = 1/(-1j*k0*Z_0)*(curlE_inc_phi+curlE_scatt_phi); % along z 
end