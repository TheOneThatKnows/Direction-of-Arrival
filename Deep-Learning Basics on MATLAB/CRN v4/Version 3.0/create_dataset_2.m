%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

custom_layer_path = ['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\' ...
    'Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers'];
addpath(custom_layer_path);

load CRN_Network_v3_0.mat
CRN_NET_1 = net;
clear net

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

%% Dataset Preparation CRN2 (2 Source)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = 1200000;

delta_phi = 1;
phi_min = 30;
phi_max = 150;
Q = (phi_max - phi_min) / delta_phi + 1;

features = zeros(2 * N - 1, numOfData);
labels = zeros(Q, numOfData);

for idx = 1:numOfData
    K = 2;
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
    labels(inds, idx) = 100;

    A_ohm = DOA.Array_Manifold(0.5, sensor_locations, doa);
    L = 70; % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');
    normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

    input1(:, :, 1) = real(normalized_R_ohm);
    input1(:, :, 2) = imag(normalized_R_ohm);
    input1(:, :, 3) = angle(R_ohm) / pi;

    features(:, idx) = predict(CRN_NET_1, input1).';
end

%% Remove Version 2 Path

rmpath(custom_layer_path);