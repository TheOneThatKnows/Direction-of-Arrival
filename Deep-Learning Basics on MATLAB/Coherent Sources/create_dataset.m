%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\');
DOA = FunctionsOfDOA();

%% Sensor Properties

% sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors
sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

%% Estimating the Number of Total Sources

numOfData = 1200000;
feature = zeros(M, M, 2);
features = zeros(M, M, 2, numOfData);

L = 70;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
labels = zeros(1, numOfData);

R_ = zeros(M, M, (M+1)/2);
R_m = zeros(M);

C_ = zeros(M, M, (M+1)/2);
C_m = zeros(M);

for idx = 1:numOfData
    K = randi(N-1);     % # of sources
    K_coherent = randi(K);

    labels(idx) = K;

    K_coherent = K_coherent * sign(K_coherent - 1);

    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

    A = DOA.Array_Manifold(sensor_locations, doa);
    s = [DOA.Coherent_Source_Generate(K_coherent, L); 
        DOA.Source_Generate(K - K_coherent, L)];
    shuffledIndices = randperm(K);
    s = s(shuffledIndices, :);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    for m = 0:(M+1)/2
    J = zeros(M - m);
    for j = 1:M-m
        J(j, M-m-j+1) = 1;
    end

    % R
    r_m = [zeros(1, m) (R(m+1:M, m+1).' * J) R(m+1, m+2:M) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        R_m(i+1, :) = r_m(start_idx:end_idx).';
    end

    R_(:, :, m+1) = R_m;

    % C
    c_m = [zeros(1, m) ones(1, 2*M-2*m-1) zeros(1, m)].';

    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        C_m(i+1, :) = c_m(start_idx:end_idx).';
    end

    C_(:, :, m+1) = C_m;
    end

    R_sum = zeros(M);
    C_sum = zeros(M);
    for m = 0:(M-1)/2
        R_sum = R_sum + R_(:, :, m+1);
        C_sum = C_sum + C_(:, :, m+1);
    end

    R_out = R_sum ./ C_sum;

    R_norm = (R_out - mean(R_out(:))) / std(R_out(:));

    feature(:, :, 1) = real(R_norm);
    feature(:, :, 2) = imag(R_norm);

    features(:, :, :, idx) = feature;
end
labels = categorical(labels);