%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load First_Level.mat

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

%% Dataset Preparation CRN2 (# of Sources = K)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = 1200000;

delta_phi = 1;
phi_min = 30;
phi_max = 150;
Q = (phi_max - phi_min) / delta_phi + 1;

features = zeros((M-1) * (M-1), numOfData);
feature = zeros(M * M, 1);

W1 = first_level.Layers(2, 1).Weights;
B1 = first_level.Layers(2, 1).Bias;

K = 2;
for idx = 1:numOfData
    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

    A_ohm = DOA.Array_Manifold(sensor_locations, doa);
    L = 50; % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');

    real_R = real(R_ohm);
    imag_R = imag(R_ohm);
    
    for i = 1:M
        feature(i) = real_R(i, i);
    end
    ind = M + 1;
    for i = 2:M
        for j = 1:i-1
            feature(ind) = real_R(i, j);
            feature(ind + 1) = imag_R(i, j);

            ind = ind + 2;
        end
    end
    zero_mean_feature = feature - mean(feature);
    feature = zero_mean_feature / max(abs(zero_mean_feature));

    features(:, idx) = W1 * feature + B1;
end
labels = features;