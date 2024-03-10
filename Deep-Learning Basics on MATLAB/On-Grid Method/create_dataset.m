%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Sensor Properties

sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(M) + 1;

%% Dataset Preparation I

numOfData = 3100000;
features = zeros(2 * M, numOfData);

K = 2;          % # of sources
L = 1;          % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
labels = zeros((phi_max - phi_min)/delta_phi + 1, numOfData);
for idx = 1:numOfData
    doa = (-2 * delta_phi) * ones(1, K);
    i = 1;
    while true
        temp_angle = phi_min + rand * (phi_max - phi_min);
        temp_array = abs(doa - temp_angle);
        if any(temp_array < delta_phi)
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
    doa = round(doa);

    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    idx_label = doa - phi_min + 1;
    labels(idx_label, idx) = abs(s);
    SNR_dB = 30 * rand;     % min = 0 dB, max = 30 dB
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    features(:, idx) = [real(y); imag(y)];
end