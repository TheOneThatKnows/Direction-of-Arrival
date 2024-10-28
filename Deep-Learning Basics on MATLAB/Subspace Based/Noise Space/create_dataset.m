%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

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

features = zeros(2 * M, numOfData);
labels = zeros(Q, numOfData);

K = 2;
for idx = 1:numOfData
    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
    rounded_doa = round(doa);

    inds = (rounded_doa - phi_min) / delta_phi + 1;
    labels(inds, idx) = 10;

    A_ohm = DOA.Array_Manifold(sensor_locations, doa);
    L = 50; % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');
    
    [eig_vecs, ~] = eig(R_ohm);
    noise_vecs = eig_vecs(:, 1:M-K);
    G = noise_vecs * noise_vecs';

    G_u = [real(G(:, 1)); imag(G(:, 1))];
    mean_G = mean(G_u);

    features(:, idx) = (G_u - mean_G) / max(abs(G_u - mean_G));
end