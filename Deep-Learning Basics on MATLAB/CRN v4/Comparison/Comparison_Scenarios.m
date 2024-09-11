%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival' ...
    '\Deep-Learning Basics on MATLAB\Custom Layers']);

run Load_Networks.m % Load The Networks

NET = Network_Functions();

%% 
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

SNR_dB_vals = -10:1:10;
EPOCHS = 5000;

noOfMethods = 14;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angles = phi_min:delta_phi:phi_max;

% Below is used for Sparse 1D Network
A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angles);

for epoch = 1:EPOCHS
    while true
        doa = phi_min + rand(1, 2) * (phi_max - phi_min);
        if abs(doa(2) - doa(1)) > delta_phi
            break
        end
    end
    doa = sort(doa);

    s = DOA.Source_Generate(K, L);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y = A * s + n;

        % CBF
        spec = CBF(angles, sensor_locations, y);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(1, idx) = RMSE(1, idx) + rmse(doa_est, doa);

        Ry = (1 / L) * (y * y');

        % Capon
        spec = Capon(angles, sensor_locations, y, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(2, idx) = RMSE(2, idx) + rmse(doa_est, doa);

        % MUSIC
        spec = MUSIC(angles, sensor_locations, Ry, M, K);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(3, idx) = RMSE(3, idx) + rmse(doa_est, doa);

        % CRN_2 v1
        R = NET.CRN2_Function_v1(CRN_v1, M, Ry);
        spec = MUSIC(angles, 0:N-1, R, N, K);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(4, idx) = RMSE(4, idx) + rmse(doa_est, doa);

        % CRN_2 v2.1
        spec = NET.CRN2_Function_v2_1(CRN_v2, CRN_v2_1, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(5, idx) = RMSE(5, idx) + rmse(doa_est, doa);

        % CRN_2 v2.2
        spec = NET.CRN2_Function_v2_2(CRN_v2_2, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(6, idx) = RMSE(6, idx) + rmse(doa_est, doa);

        % CRN_2 v2.4
        spec = NET.CRN2_Function_v2_4(CRN_v2_4, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(7, idx) = RMSE(7, idx) + rmse(doa_est, doa);

        % CRN_2 v2.5
        spec = NET.CRN2_Function_v2_5(CRN_v2_5, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(8, idx) = RMSE(8, idx) + rmse(doa_est, doa);

        % CRN_2 v2.6
        spec = NET.CRN2_Function_v2_6(CRN_v2_6, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(9, idx) = RMSE(9, idx) + rmse(doa_est, doa);

        % CRN_2 v2.7
        spec = NET.CRN2_Function_v2_7(CRN_v2_7, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(10, idx) = RMSE(10, idx) + rmse(doa_est, doa);

        % CRN_2 v2.8
        spec = NET.CRN2_Function_v2_8(CRN_v2_8, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(11, idx) = RMSE(11, idx) + rmse(doa_est, doa);

        % DNN v2
        spec = NET.SparseDNN_Function_v2(DNN_v2, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(12, idx) = RMSE(12, idx) + rmse(doa_est, doa);

        % DNN v3
        spec = NET.SparseDNN_Function_v3(DNN_v3, N, Ry, DOA, sensor_locations);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(13, idx) = RMSE(13, idx) + rmse(doa_est, doa);

        % Sparse 1D
        spec = NET.Sparse_1D_Function(Sparse_1D_Net, Ry, DOA, A_sparse);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(14, idx) = RMSE(14, idx) + rmse(doa_est, doa);
    end
    if rem(epoch, 10) == 0
        disp(epoch)
    end
end

% RMSE = (1 / EPOCHS) * RMSE;
temp_RMSE = RMSE;
RMSE = (1 / (epoch - 1)) * RMSE;

%% Plot 

clear; clc; close all;
load RMSE_vals.mat
SNR_dB_vals = -10:1:10;

figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(14, :));
plot(SNR_dB_vals, RMSE(12, :));
plot(SNR_dB_vals, RMSE(13, :));
plot(SNR_dB_vals, RMSE(4, :));
plot(SNR_dB_vals, RMSE(5, :));
plot(SNR_dB_vals, RMSE(6, :));
plot(SNR_dB_vals, RMSE(7, :));
plot(SNR_dB_vals, RMSE(8, :));
plot(SNR_dB_vals, RMSE(9, :), '-*');
plot(SNR_dB_vals, RMSE(10, :));
plot(SNR_dB_vals, RMSE(11, :));

xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', '1D-Sparse Net', 'DNN v2', 'DNN v3', ...
    'CRN2 v1', 'CRN2 v2.1', 'CRN2 v2.2' , 'CRN2 v2.4', 'CRN2 v2.5', ...
    'CRN2 v2.6', 'CRN v2.7', 'CRN v2.8');
title('SNR vs RMSE')

%% Functions

% DOA Estimator
function doa_est = DOA_Estimator(spec, angles)
spec = [0 spec 0];
[mags, inds] = findpeaks(spec);
doa_est = zeros(1, 2);
[~, ind] = max(mags);
idx = inds(ind);
doa_est(1) = angles(idx - 1);
mags = [mags(1:ind-1) mags(ind+1:end)];
inds = [inds(1:ind-1) inds(ind+1:end)];
[~, ind] = max(mags);
idx = inds(ind);
if isempty(idx)
    doa_est(2) = doa_est(1);
else
    doa_est(2) = angles(idx - 1);
end

doa_est = sort(doa_est);
end

% CBF
function spec = CBF(angles, sensor_locations, y)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = abs(h' * (y * y') * h);
end
spec = 10 * log10(spec);
end

% Capon
function spec = Capon(angles, sensor_locations, y, Ry)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    h = (Ry \ a) / (a' / Ry * a);
    spec(i) = abs(h' * (y * y') * h);
end
end

% MUSIC
function spec = MUSIC(angles, sensor_locations, Ry, M, n)
[eig_vecs, ~] = eig(Ry);
U_N = eig_vecs(:, 1:M-n);
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = 1 / abs(a' * (U_N * U_N') * a);
end
end