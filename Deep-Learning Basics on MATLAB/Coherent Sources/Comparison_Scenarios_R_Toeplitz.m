%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

rng('shuffle')

%% R_Toeplitz, R_Toeplitz2, R_Toeplitz3 have the same accuracy 

sensor_locations = 0:10; % ULA with 11 sensors
M = length(sensor_locations);
K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

SNR_dB_vals = -10:1:20;
EPOCHS = 5000;

noOfMethods = 4; % R_Toeplitz, R_Toeplitz2, R_Toeplitz3, R_Toeplitz4
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angle_spec = phi_min:delta_phi/10:phi_max;

for epoch = 1:EPOCHS
    doa = DOA.DOA_Generate(K, phi_min, phi_max, 2 * delta_phi);
    A = DOA.Array_Manifold(sensor_locations, doa);

    vars = ones(K, 1);
    vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
    s = DOA.Source_Generate_Final(K, K_coherent, L, vars);

    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y = A * s + n;
        Ry = (1 / L) * (y * y');

        method = 1;

        % R_Toeplitz + MUSIC
        R_toeplitz = R_Toeplitz(Ry, "full");
        spec = DOA.MUSIC(K, R_toeplitz, sensor_locations, angle_spec);
        doa_est = DOA.DOA_Estimate(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % R_Toeplitz2 + MUSIC
        R_toeplitz = R_Toeplitz2(Ry, "full");
        spec = DOA.MUSIC(K, R_toeplitz, sensor_locations, angle_spec);
        doa_est = DOA.DOA_Estimate(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % R_Toeplitz3 + MUSIC
        R_toeplitz = R_Toeplitz3(Ry, "full");
        spec = DOA.MUSIC(K, R_toeplitz, sensor_locations, angle_spec);
        doa_est = DOA.DOA_Estimate(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % R_Toeplitz4 + MUSIC
        R_toeplitz = R_Toeplitz4(Ry, "full");
        spec = DOA.MUSIC(K, R_toeplitz, sensor_locations, angle_spec);
        doa_est = DOA.DOA_Estimate(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
    end
    if rem(epoch, 10) == 0
        disp(epoch + " / " + EPOCHS)
    end
end

RMSE = (1 / EPOCHS) * RMSE;
% RMSE = (1 / epoch) * RMSE;

figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(4, :));
% plot(SNR_dB_vals, 10*log10(RMSE(1, :)), 'b--o');
% plot(SNR_dB_vals, 10*log10(RMSE(2, :)), 'r*');
% plot(SNR_dB_vals, 10*log10(RMSE(3, :)));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('Toeplitz ULA MUSIC', 'Toeplitz Sparse MUSIC', 'SS-MUSIC');
title_text = "SNR vs RMSE (K=" + K + ", K_{coherent}="+ K_coherent + ")";
title(title_text)

ymax = ceil(max(RMSE(:)));
ylim([0 ymax])