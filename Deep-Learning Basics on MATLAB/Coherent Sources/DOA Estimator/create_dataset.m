%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Coherent Sources']);

%% Sensor Properties

sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) - sensor_locations(1) + 1;

%% Prepare Dataset

numOfData = 1200000;
feature = zeros(M, M, 3);
features = zeros(M, M, 3, numOfData);

L = 70;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
Q = (phi_max - phi_min) / delta_phi + 1;
labels = zeros(numOfData, Q);

K = 3;
K_coherent = 2;
for idx = 1:numOfData
    % K = randi(N-1);     % # of sources
    % K_coherent = randi(K);
    % K_coherent = K_coherent * sign(K_coherent - 1);

    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
    for i = 1:K
        c = floor((doa(i) - phi_min) / delta_phi);

        phi_a = phi_min + c * delta_phi;
        phi_b = phi_a + delta_phi;

        percentages = [phi_a phi_b; 1 1] \ [doa(i)*100; 100];

        idx1 = c + 1;
        idx2 = c + 2;
        labels(idx, idx1:idx2) = labels(idx, idx1:idx2) + percentages.';
    end

    A = DOA.Array_Manifold(sensor_locations, doa);
    s = [DOA.Coherent_Source_Generate(K_coherent, L); 
        DOA.Source_Generate(K - K_coherent, L)];
    shuffledIndices = randperm(K);
    s = s(shuffledIndices, :);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    R_toeplitz = R_Toeplitz(R, "half");

    R_norm = (R_toeplitz - mean(R_toeplitz(:))) / std(R_toeplitz(:));

    feature(:, :, 1) = real(R_norm);
    feature(:, :, 2) = imag(R_norm);
    feature(:, :, 3) = angle(R_norm) / pi;

    features(:, :, :, idx) = feature;
end