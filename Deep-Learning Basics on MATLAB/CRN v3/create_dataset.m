%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

%% Dataset Preparation

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = 1;
phi_min = 30;
phi_max = 150;
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
    L = 100; % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    r = abs(R(:));
    max_r = max(r);

    features(:, :, 1, idx) = (real(R)) / max_r;
    features(:, :, 2, idx) = (imag(R)) / max_r;
end