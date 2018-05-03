clear all; close all; clc;

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

H_inc_0 = 1/Z_0;

% Anonymous Functions: sqr = @(x) x.^2; sqr(5) % we get 25
G0 = @(k0r) (1/4j)*besselh(0,2,k0r);
exp_gamma = 1.7810724179901979852; % == exp(0.577215664901532 the Euler�Mascheroni constant)

%% plot |G0|
% k0r = linspace(0,4*pi, 500);
% G0_k0r = G0(k0r);
% figure(1); 
% hold on;
% plot(k0r,abs (G0_k0r), 'DisplayName', '|G0|');
% plot(k0r,real(G0_k0r), 'DisplayName', 'real(G0)');
% plot(k0r,imag(G0_k0r), 'DisplayName', 'imag(G0)');
% legend('show'); % kind of different in 2018a

%%
r = 1*lambda_;
N = 500;
theta = linspace(0, 2*pi, N+1).'; theta(end) = []; % colum
d_theta = 2*pi/N;
s_n = r*d_theta;
% [-d_theta/2, +d_theta/2], ...
xx = r*cos(theta); yy = r*sin(theta);

R = sqrt((xx - xx.').^2 + (yy - yy.').^2); % |rho_m - rho_n|
norm_n_x = cos(theta).';
norm_n_y = sin(theta).';
term2 = ((xx - xx.').*norm_n_x + (yy - yy.').*norm_n_y) ./ R;  % m as | ; n as --
Zmn = (k0*s_n/4j)*besselh(1,2,k0*R).*term2;

diag_value = -1/2;
Zmn(1:1+N:end) = diag_value;

Vm = 1/Z_0*exp(-1j*k0*xx); % add a phase shift

J = Zmn \ Vm;
figure(2); hold on;
plot(theta*(180/pi),abs(J)/H_inc_0); xlim([0, 360]); xlabel('\phi (degree)');ylabel('J_s / H_0')
xticks([0, 90, 180, 270, 360]);


%% Hz else where
N2 = 1000+1;
x_out = lambda_*linspace(-5,+5,N2).';
y_out = lambda_*linspace(-5,+5,N2).';
H_out = complex(zeros(N2, N2));
h = lambda_*(+5-(-5))/N2;

for jj = 1:N2
    y_out0 = y_out(jj);
    mask = sqrt(x_out.^2 + y_out0.^2) <= r+0.5*h;
 
    R = sqrt((x_out-xx.').^2 + (y_out0-yy.').^2);
    % full field
    term2 = ((x_out - xx.').*norm_n_x + (y_out0 - yy.').*norm_n_y) ./ R;
    % H_out(:, jj)= 1/Z_0*exp(-1j*k0*x_out) + (besselh(1,2,k0*R).*term2)*J*(-k0/4j*s_n);
    % scattered field
    H_out(:, jj)= (besselh(1,2,k0*R).*term2)*J*(-k0/4j*s_n);
    aaa=besselh(1,2,k0*R);
    H_out(mask, jj) = NaN;
end

figure(3); hold on;
surf(x_out, y_out, abs(H_out.')*Z_0, 'EdgeColor','none','LineStyle','none','FaceLighting','gouraud'); 
% The 'phong' value has been removed. Use 'gouraud' instead.
view(2);xlim([-5, +5]);ylim([-5, +5]);xticks([-5:1:5]);yticks([-5:1:5]);colorbar;axis('image');



