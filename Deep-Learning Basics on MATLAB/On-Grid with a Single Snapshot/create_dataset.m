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

ula_2_sparse = [1 0 0 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0 0;
                0 0 0 0 0 0 0 1 0 0;
                0 0 0 0 0 0 0 0 0 1];

numOfData = 1200000;

delta_phi = 1;
phi_min = 30;
phi_max = 150;
Q = (phi_max - phi_min) / delta_phi + 1;

features = zeros(2 * M, numOfData);
labels = zeros(2 * N, numOfData);

K = 2; 
L = 1; % # of snapshots
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

    A = DOA.Array_Manifold(0:N-1, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s;
    y_ohm = ula_2_sparse * y + n;

    features(:, idx) = [real(y_ohm); imag(y_ohm)];
    labels(:, idx) = [real(y); imag(y)];
end