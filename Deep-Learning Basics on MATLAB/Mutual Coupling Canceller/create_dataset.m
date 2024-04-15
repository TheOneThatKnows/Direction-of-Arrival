%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [1 2 5 8 10] - 1; % SLA with 5 sensors

%% Dataset Preparation CRN2 (2 Source)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = 3100000;
features = zeros(M*M, numOfData);
labels = zeros(2*M-1, numOfData);

K = 2;
L = 100;
delta_phi = 1;
phi_min = 30;
phi_max = 150;
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

    C = DOA.Mutual_Coupling(N, 0.5, M, sensor_locations);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = 0 + rand * 30;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = C * A * s + n;
    Ry = (1 / L) * (y * y');

    Rs = (1 / L) * (s * s');
    R = A * Rs * A';

    labels(1:M, idx) = real(R(:, 1));
    labels(M+1:end, idx) = imag(R(2:end, 1));

    re_R = real(Ry);
    im_R = imag(Ry);

    features(1:M, idx) = diag(re_R);
    ind = M + 1;
    for i = 2:M
        for j = 1:i-1
            features(ind:ind+1, idx) = [re_R(i, j); im_R(i, j)];
            ind = ind + 2;
        end
    end
end