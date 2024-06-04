%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = 0:7; % ULA with 8 sensors
% sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors
% sensor_locations = [0 1 2 6 10 13]; % SLA with 6 sensors

%% Dataset Preparation

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = 1;
phi_min = 60;
phi_max = 119;
angles = phi_min:delta_phi:phi_max;
Q = length(angles);

numOfData = 15800;
% numOfData = 200;
features = zeros(M, M, 2, numOfData);
labels = zeros(Q, numOfData);

for idx = 1:numOfData
    doa = (-2 * delta_phi) * ones(1, K);
    i = 1;
    while true
        temp_angle = phi_min + rand * (phi_max - phi_min);
        temp_array = abs(doa - temp_angle);
        if any(temp_array < delta_phi * 2)
            i = i - 1;
        else
            doa(i) = temp_angle;
        end
        if i == K
            break
        end
        i = i + 1;
    end
    doa = sort(doa);
    rounded_doa = round(doa);
    inds = (rounded_doa - phi_min) / delta_phi + 1;
    labels(inds, idx) = 1;

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    L = 256; % # of snapshots
    Tau = randi(20); % # delay
    s_temp = DOA.Source_Generate(1, L + Tau);
    s = zeros(K, L);
    s(1, :) = s_temp(1:L);
    s(2, :) = s_temp(Tau+1:end);
    SNR_dB = 0 - rand * 20;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    r = abs(R(:));
    % mean_r = mean(r);
    % var_r = var(r);
    max_r = max(r);

    % features(:, :, 1, idx) = (real(R) - mean_r) / var_r;
    % features(:, :, 2, idx) = (imag(R) - mean_r) / var_r;

    features(:, :, 1, idx) = (real(R)) / max_r;
    features(:, :, 2, idx) = (imag(R)) / max_r;
end