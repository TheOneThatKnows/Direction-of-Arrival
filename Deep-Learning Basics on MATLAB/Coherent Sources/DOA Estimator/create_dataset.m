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

% sensor_locations = [1 2 3 4 9 13] - 1; % NA v2 with 6 sensors
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

%% Estimating the Number of Total Sources

numOfData = 1200000;
feature = zeros(N, N, 2);
features = zeros(N, N, 2, numOfData);

L = 70;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
Q = (phi_max - phi_min) / delta_phi + 1;
labels = zeros(numOfData, Q);

K = 5;
K_coherent = 2;
for idx = 1:numOfData
    % K = randi(N-1);     % # of sources
    % K_coherent = randi(K);
    % K_coherent = K_coherent * sign(K_coherent - 1);

    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
    doa_inds = (round(doa) - phi_min) / delta_phi + 1;
    labels(idx, doa_inds) = 100;

    A = DOA.Array_Manifold(sensor_locations, doa);
    s = [DOA.Coherent_Source_Generate(K_coherent, L); 
        DOA.Source_Generate(K - K_coherent, L)];
    shuffledIndices = randperm(K);
    s = s(shuffledIndices, :);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    % R_out = Spatial_Smoothing(DOA, sensor_locations, R);

    z = R(:);
    z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
    R_out = zeros(N);
    for i = 1:N
        z1_i = z1(i:i + N - 1);
        R_out = R_out + (1 / N) * (z1_i * z1_i');
    end

    R_norm = (R_out - mean(R_out(:))) / std(R_out(:));

    feature(:, :, 1) = real(R_norm);
    feature(:, :, 2) = imag(R_norm);

    features(:, :, :, idx) = feature;
end