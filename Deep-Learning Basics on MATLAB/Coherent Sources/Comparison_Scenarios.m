%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Coherent Sources']);

rng('shuffle')

%% 
sensor_locations_ula = 0:4; % ULA with 5 sensors
sensor_locations_mra = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations_mra);
N = sensor_locations_mra(M) + 1;
K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

SNR_dB_vals = -10:1:10;
EPOCHS = 5000;

noOfMethods = 3; % R_Toeplitz+MUSIC for ULA, R_Sparse_Toeplitz+MUSIC for MRA, SS-MUSIC for MRA, R_Sparse_Toeplitz+SS-MUSIC for MRA
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angle_spec = phi_min:delta_phi/10:phi_max;

for epoch = 1:EPOCHS
    doa = DOA.DOA_Generate(K, phi_min, phi_max, 2 * delta_phi);

    s = [DOA.Coherent_Source_Generate(K_coherent, L);
        DOA.Source_Generate(K - K_coherent, L)];
    shuffledIndices = randperm(K);
    s = s(shuffledIndices, :);
    A_ula = DOA.Array_Manifold(sensor_locations_ula, doa);
    A_mra = DOA.Array_Manifold(sensor_locations_mra, doa);
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y_ula = A_ula * s + n;
        y_mra = A_mra * s + n;

        method = 1;

        % R_Toeplitz + MUSIC
        Ry_ula = (1 / L) * (y_ula * y_ula');
        R_toeplitz = R_Toeplitz(Ry_ula, "half");
        spec = DOA.MUSIC(K, R_toeplitz, sensor_locations_ula, angle_spec);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;


        Ry_mra = (1 / L) * (y_mra * y_mra');

        % Sparse-Toeplitz + MUSIC
        R_sparse_toeplitz = Sparse_Toeplitz(DOA, Ry_mra, sensor_locations_mra);
        spec = DOA.MUSIC(K, R_sparse_toeplitz, sensor_locations_mra, angle_spec);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % SS-MUSIC
        spec = DOA.SS_MUSIC(K, Ry_mra, sensor_locations_mra, angle_spec);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        doa_est = sort(doa_est);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
    end
    if rem(epoch, 10) == 0
        disp(epoch + " / " + EPOCHS)
    end
end

RMSE = (1 / EPOCHS) * RMSE;
% RMSE = (1 / epoch) * RMSE;

%%
figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
% plot(SNR_dB_vals, 10*log10(RMSE(1, :)), 'b--o');
% plot(SNR_dB_vals, 10*log10(RMSE(2, :)), 'r*');
% plot(SNR_dB_vals, 10*log10(RMSE(3, :)));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('Toeplitz ULA MUSIC', 'Toeplitz Sparse MUSIC', 'SS-MUSIC');
title_text = "SNR vs RMSE (K=" + K + ", K_{coherent}="+ K_coherent + ")";
title(title_text)
% ylim([0 25])
%% Functions

% DOA Estimator

function doa_est = DOA_Estimator(spec, angles, K)
spec = [-inf spec -inf];
[mags, mags_inds] = findpeaks(spec);
doa_est = zeros(1, K);
[~, sorted_inds] = sort(mags, "descend");
mags_inds = mags_inds(sorted_inds);
mags_inds = [mags_inds zeros(1, K-1)];
idx = mags_inds(1);
doa_est(1) = angles(idx - 1);

for i = 2:K
    idx = mags_inds(i);
    if idx == 0
        doa_est(i:K) = doa_est(i-1);
        break
    else
        doa_est(i) = angles(idx - 1);
    end
end

doa_est = sort(doa_est);
end

% DOA Estimator Off-Grid
function doa_est = DOA_Estimator_Off_Grid(spec, angles)
spec = [0 spec 0];
angles = [0 angles 0];
[mags, mags_inds] = findpeaks(spec);
doa_est = zeros(1, 3);
[~, sorted_inds] = sort(mags, "descend");
mags_inds = mags_inds(sorted_inds);
mags_inds = [mags_inds zeros(1, 2)];

idx = mags_inds(1);
doa_est(1) = ([angles(idx-1) angles(idx) angles(idx+1)] * [spec(idx-1) spec(idx) spec(idx+1)].') / sum(spec(idx-1:idx+1));

idx = mags_inds(2);
if idx == 0
    doa_est(2:3) = doa_est(1);
    return
end
doa_est(2) = ([angles(idx-1) angles(idx) angles(idx+1)] * [spec(idx-1) spec(idx) spec(idx+1)].') / sum(spec(idx-1:idx+1));

idx = mags_inds(3);
if idx == 0
    c = randi(2);
    doa_est(3) = (c - 1) * doa_est(1) +  rem(c, 2) * doa_est(2);
    doa_est = sort(doa_est);
    return
end
doa_est(3) = ([angles(idx-1) angles(idx) angles(idx+1)] * [spec(idx-1) spec(idx) spec(idx+1)].') / sum(spec(idx-1:idx+1));

doa_est = sort(doa_est);
end

% DL-Network

function spec = Net_Function(net, M, R)
R_norm = (R - mean(R(:))) / std(R(:));

feature = zeros(M, M, 3);
feature(:, :, 1) = real(R_norm);
feature(:, :, 2) = imag(R_norm);
feature(:, :, 3) = angle(R_norm) / pi;

spec = predict(net, feature).';
spec = spec / max(spec);
end