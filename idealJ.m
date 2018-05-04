% D(hankelh2[n,x],x)
% 
% D(BesselJ[n,x],x)

if strcmp(mode, 'TE Hz')
bn = @(n, ka) -(1j)^(-n)*( ...
  ( 0.5*(besselj(n-1,ka)   - besselj(n+1,ka))   )/ ...
  ( 0.5*(besselh(n-1,2,ka) - besselh(n+1,2,ka)) )  ...
  );
ka = 2*pi;

maxN = 40;
b_n = zeros(2*maxN+1,1);
for n=-maxN:maxN
    b_n(n+maxN+1) = bn(n, ka);
end




coeff = zeros(2*maxN+1, 1);
for n=-maxN:maxN
    coeff(n+maxN+1) = -b_n(n+maxN+1)*besselh(n,2,ka);
end
exp_theta = exp(1j* (-maxN:maxN).*theta);

J = exp_theta * coeff - exp(-1j*2*pi*cos(theta));

figure(2);
plot(theta, abs(J));
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
    xlabel('index n');title('a_n \cdot H^{(2)}_n(k\rho) | \rho = radius');
    
    exp_theta = exp(1j* (-maxN:maxN).*theta);
    curlE_scatt_phi = -exp_theta*(a_n_H2p*k0); 
    % - sign due to 
    % \nabla \times E^{scatt} = \frac{1}{\rho} \diff{E_z}{\phi} \hat{\rho}
    % + -\diff{E_z}{\rho} \hat{\phi}
    J_ideal = 1/(-1j*k0*Z_0)*(curlE_inc_phi+curlE_scatt_phi); % along z 
end