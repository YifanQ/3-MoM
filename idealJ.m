% D(hankelh2[n,x],x)
% 
% D(BesselJ[n,x],x)

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

figure(1);
plot(-maxN:maxN, real(b_n), -maxN:maxN, imag(b_n));

N = 500;
theta = linspace(0, 2*pi, N+1).'; theta(end) = []; % colum

coeff = zeros(2*maxN+1, 1);
for n=-maxN:maxN
    coeff(n+maxN+1) = -b_n(n+maxN+1)*besselh(n,2,ka);
end
exp_theta = exp(1j* (-maxN:maxN).*theta);

J = exp_theta * coeff - exp(-1j*2*pi*cos(theta));

figure(2);
plot(theta, abs(J));