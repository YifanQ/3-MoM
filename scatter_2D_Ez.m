clear all; close all; clc;

plotE_out = false;

%% EM parameters ==================================================
epsilon_    =8.854e-012; %8.854 187 817
mu_         =1.257e-006; %1.256 637 061
%Speed of light
c_=1/sqrt(epsilon_*mu_); %299 792 458
%Impedance of free space
Z_0=sqrt(mu_/epsilon_); %Z_0 = E/H = 376.730 313 Ohm

lambda_     = 1.0;
frequency   = c_/lambda_;

omega       = 2*pi*frequency;
k0          = 2*pi/lambda_; % = 2*pi*frequency/c_ = omega/c_;

E_inc_0 = 1.0;
H_inc_0 = E_inc_0/Z_0;

% Anonymous Functions: sqr = @(x) x.^2; sqr(5) % we get 25
G0 = @(k0r) (1/4j)*besselh(0,2,k0r); %-% Equ 10.2.4
exp_gamma = 1.7810724179901979852; % == exp(0.577215664901532 the Euler–Mascheroni constant)

%% plot |G0|
% k0r = linspace(0,4*pi, 500);
% G0_k0r = G0(k0r);
% figure(1);
% hold on;
% plot(k0r,abs (G0_k0r), 'DisplayName', '|G0|');
% plot(k0r,real(G0_k0r), 'DisplayName', 'real(G0)');
% plot(k0r,imag(G0_k0r), 'DisplayName', 'imag(G0)');
% legend('show'); % kind of different in 2018a

%% define geometry
N = 500;
theta = linspace(0, 2*pi, N+1).'; theta(end) = []; % column
d_theta = 2*pi/N;

r = 1*lambda_;
s_n = r*d_theta;
% [-d_theta/2, +d_theta/2], ...
xx = r*cos(theta); yy = r*sin(theta);

R = sqrt((xx - xx.').^2 + (yy - yy.').^2); % |rho_m - rho_n|

%-% Equ 10.2.25
Zmn = (k0*Z_0*s_n/4)*besselh(0,2,k0*R);
diag_value = (k0*Z_0*s_n/4)*(  1-1j*(2/pi)*log(k0*s_n*exp_gamma/( 4*exp(1.0) ))  );
Zmn(1:1+N:end) = diag_value;

Vm = E_inc_0*exp(-1j*k0*xx);

J = Zmn \ Vm;
mode = 'TM Ez';
idealJ;

figure(2);
ax1=subplot(4,1,[1 2 3]); hold on;
plot(theta*(180/pi),abs(J)/H_inc_0); xlim([0, 360]); xlabel('\phi (degree)');ylabel('J_s / H_0')
xticks([0, 90, 180, 270, 360]);

figure(2);
ax2=subplot(4,1,4); hold on;
plot(theta*(180/pi),( abs(J-J_ideal) ) / H_inc_0); xlim([0, 360]); xlabel('\phi (degree)');ylabel('|J_s-J^*_s| / H_0')
xticks([0, 90, 180, 270, 360]);

J_error = norm(J-J_ideal)/norm(J_ideal);
subplot(4,1,[1 2 3]);
title(sprintf('surface J_z, # of unknowns = %d, error ||J-J^*||_2 / ||J^*||_2 = %0.3f', N, J_error))

set([ax1 ax2], 'FontSize', 16);

%% Ez at rho for sigma_2D
rho_list = lambda_*[10, 100, 1000];
figure(10); hold on;
for ii = 1:length(rho_list)
    rho = rho_list(ii);

    N1 = 1000;
    theta1 = linspace(0, 2*pi, N1+1).'; theta1(end) = []; % column

    x_out = rho*cos(theta1); y_out = rho*sin(theta1);

    R = sqrt((x_out-xx.').^2 + (y_out-yy.').^2);
    % scattered field
    E_scatt = G0(k0*R)*J*(-1j*k0*Z_0*s_n);
    sigma_2D = 2*pi*rho*abs( E_scatt ).^2 / abs(E_inc_0).^2;
    % plot(theta1*(180/pi),sigma_2D, 'DisplayName', sprintf('\\rho = %d',rho));
    plot(theta1*(180/pi),10*log10(sigma_2D/lambda_), 'DisplayName', sprintf('\\rho = %d \\lambda',rho));
end
xlim([0, 360]); xlabel('\phi (degree)');
% ylabel('\sigma_{2D}');
ylabel('\sigma_{2D}/\lambda (dB)');
xticks([0, 90, 180, 270, 360]);
legend('show');

if ~plotE_out
    return
end

%% Ez else where
N2 = 1000+1;
x_out = lambda_*linspace(-5,+5,N2).';
y_out = lambda_*linspace(-5,+5,N2).';
E_scatt = complex(zeros(N2, N2));
E_full = complex(zeros(N2, N2));

for jj = 1:N2
    y_out0 = y_out(jj);
    mask = sqrt(x_out.^2 + y_out0.^2) <= r;

    R = sqrt((x_out-xx.').^2 + (y_out0-yy.').^2);
    % scattered field
    E_scatt(:, jj) = G0(k0*R)*J*(-1j*k0*Z_0*s_n);
    % full field
    E_full(:, jj)  = E_scatt(:, jj) + E_inc_0*exp(-1j*k0*x_out); % E_xxxx(x_idx, y_idx)

    E_scatt(mask, jj) = NaN;
    E_full (mask, jj) = NaN;
end

plotMat(3, x_out, y_out, abs(E_scatt.')/E_inc_0);
plotMat(4, x_out, y_out, abs(E_full.' )/E_inc_0);
plotMat(5, x_out, y_out,real(E_scatt.')/E_inc_0);
plotMat(6, x_out, y_out,real(E_full.' )/E_inc_0);

function plotMat(figID, x_out, y_out, mat)
figure(figID); hold on;
surf(x_out, y_out, mat, 'EdgeColor','none','LineStyle','none','FaceLighting','gouraud');
% The 'phong' value has been removed. Use 'gouraud' instead.
view(2);xlim([-5, +5]);ylim([-5, +5]);xticks([-5:1:5]);yticks([-5:1:5]);colorbar;axis('image');
end
