%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Coherent Sources']);

rng shuffle

%% Define The Sensor Array

M = 4; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations_ca = sort([sensor_locations_1 sensor_locations_2(2:end)]);
sensor_locations_oca = [0 max(M,N)+sensor_locations_ca(2:end)];

M = length(sensor_locations_oca);
N = sensor_locations_oca(M) + 1;

% virtual sensors

virtual_sensor_locations = 0:N-1;
virtual_sensor_locations([15 17 19 20]) = 0;
virtual_sensor_locations = sort(virtual_sensor_locations, "ascend");
virtual_sensor_locations = virtual_sensor_locations(5:end);

%% Prepare Dataset

numOfData = 1200000;
feature = zeros(N, N, 2);
features = zeros(N, N, 2, numOfData);

L = 70;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
Q = (phi_max - phi_min) / delta_phi + 1;
labels = zeros(numOfData, 2*N-1);

K = 5;
for idx = 1:numOfData
    K_coherent = randi(K);
    K_coherent = K_coherent * sign(K_coherent - 1);

    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

    % Array Manifold
    A = DOA.Array_Manifold(sensor_locations_oca, doa);

    % signal generate
    vars = ones(K-(K_coherent-1)*sign(K_coherent), 1);
    s = DOA.Source_Generate_Final(K, K_coherent, L, vars);

    % noise generate
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);

    % measurements
    y = A * s + n;
    R = (1 / L) * (y * y');

    R_toeplitz = toeplitz(Virtual_Covariance_Column(DOA, R, sensor_locations_oca)');

    R_norm = (R_toeplitz - mean(R_toeplitz(:))) / std(R_toeplitz(:));

    feature(:, :, 1) = real(R_norm);
    feature(:, :, 2) = imag(R_norm);

    features(:, :, :, idx) = feature;

    A2 = DOA.Array_Manifold(virtual_sensor_locations, doa);
    v = R_toeplitz(:, 1);
    v = v(virtual_sensor_locations + 1);
    R_in = toeplitz(v');

    % incoherent source matrix
    Rs = ((A2' * A2) \ A2') * R_in * (A2 / (A2' * A2));

    A3 = DOA.Array_Manifold(0:N-1, doa);
    R_out = A3 * diag(diag(Rs)) * A3';

    labels(idx, :) = [real(R_out(:, 1).') imag(R_out(2:end, 1).')];
end