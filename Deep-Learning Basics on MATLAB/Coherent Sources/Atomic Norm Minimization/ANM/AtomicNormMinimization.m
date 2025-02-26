%% Article Name: 
% DOA Estimation of Coherent Sources Using Coprime Array via Atomic Norm Minimization 
% https://ieeexplore.ieee.org/abstract/document/9785899

%% Initialize

clear; clc; close all;
addpath('E:\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Define The Sensor Array

M = 3; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations = sort([sensor_locations_1 sensor_locations_2(2:end)]);

M = length(sensor_locations);
N = sensor_locations(M) + 1;

%% Simulate The Environment

K = 3;          % # of sources
K_coherent = 3;
L = 500;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

% doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
% doa = [50 120];
doa = [90 111.5 140];

vars = ones(K - K_coherent + 1, 1);
% vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
% alpha = [1; exp(1i * pi / 6)];
alpha = [1; exp(1i * pi / 6); exp(1i * pi / 10)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars, alpha);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

y2 = zeros(N, L);
y2(sensor_locations+1, :) = y;

Rv = (1 / L) * (y2 * y2');
J = zeros(N, M);
for i = 1:M
    J(sensor_locations(i)+1, i) = 1;
end
R_hat = Rv * J;
G = zeros(N, M);
G(sensor_locations+1, :) = ones(M);

[W, u, R0] = solve_optimization(G, R_hat);

R_final = toeplitz(u');

%% Plot The Spatial Spectrum

figure; 
ss = DOA.MUSIC(K, R_final, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss));
xlabel('angle range (deg)')
ylabel('amplitude')
title('Spatial Spectrum')

%%

L = 70;
phi_min = 30;
phi_max = 150;
delta_phi = 0.2;
Q = (phi_max - phi_min) / delta_phi + 1;

SNR_dB_vals = -10:1:10;
EPOCHS = 1000;

noOfMethods = 2; % ANM + MUSIC, FBSS + MUSIC
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angle_spec = phi_min:delta_phi:phi_max;

% ula_ca = zeros(M, N);
% for i = 1:M
%     ind = sensor_locations(i) + 1;
%     ula_ca(i, ind) = 1;
% end

J = zeros(N, M);
for i = 1:M
    J(sensor_locations(i)+1, i) = 1;
end
G = zeros(N, M);
G(sensor_locations+1, :) = ones(M);

K = 3;
K_coherent = 3;
for epoch = 1:EPOCHS
    % K = randint(2, 5)
    % K_coherent = randint(1, K+1)
    % K_coherent = K_coherent * sign(K_coherent-1)

    % DOA Angles
    % doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi)
    doa = [90 111.5 140];

    % Array Manifold
    % A = DOA.Array_Manifold(sensor_locations, doa);
    A = DOA.Array_Manifold(0:N-1, doa);

    % Signal Generate
    vars = ones(K - (K_coherent-1)*sign(K_coherent), 1);
    alpha = [1, exp(1i*pi/6), exp(1i*pi/10)].';
    s = DOA.Source_Generate_Final(K, K_coherent, L, vars, alpha);

    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        % Noise Generate
        n = DOA.Noise_Generate(SNR_dB, N, L);

        % Measurements
        y = A * s + n;
        y2 = zeros(N, L);
        y2(sensor_locations+1, :) = y(sensor_locations+1, :);
        R = (1 / L) * (y * y');
        R2 = (1 / L) * (y2 * y2');

        R_hat = R2 * J;
        [~, u, ~] = solve_optimization(G, R_hat);

        R_fb = DOA.FBSS(R, 6);

        method = 1;

        % ANM + MUSIC
        spec = DOA.MUSIC(K, toeplitz(u'), 0:N-1, angle_spec);
        doa_est = DOA.DOA_Estimate_Local(spec, angle_spec, doa, 5);
        RMSE(method, idx) = RMSE(method, idx) + sqrt(mse(doa, doa_est));
        method = method + 1;

        % FBSS + MUSIC
        spec = DOA.MUSIC(K, R_fb, 0:5, angle_spec);
        doa_est = DOA.DOA_Estimate_Local(spec, angle_spec, doa, 5);
        RMSE(method, idx) = RMSE(method, idx) + sqrt(mse(doa, doa_est));

        if rem(epoch, 10) == 0
            disp("Epoch: " + epoch + "/" + EPOCHS)
        end
    end
end
RMSE = RMSE / EPOCHS;

%% Plot RMSE

figure; hold on;
plot(SNR_dB_vals, RMSE(1, :))
plot(SNR_dB_vals, RMSE(2, :))
title("SNR vs RMSE")
xlabel('SNR (dB)')
ylabel('RMSE (degree)')
legend('ANM', 'FBSS')