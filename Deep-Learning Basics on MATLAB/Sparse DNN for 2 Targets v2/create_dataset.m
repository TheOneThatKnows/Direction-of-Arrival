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
features = zeros(M*M, numOfData);

K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
Q = (phi_max - phi_min)/delta_phi + 1;

labels = zeros(Q, numOfData);
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
    idx_label = doa - phi_min + 1;
    labels(idx_label, idx) = 100;

    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    re_R = real(R);
    im_R = imag(R);

    features(1:M, idx) = diag(re_R);
    ind = M + 1;
    for i = 2:M
        for j = 1:i-1
            features(ind:ind+1, idx) = [re_R(i, j); im_R(i, j)];
            ind = ind + 2;
        end
    end
end