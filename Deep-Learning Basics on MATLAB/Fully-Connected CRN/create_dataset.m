%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [1 2 5 8 10] - 1; % SLA with 5 sensors

%% Dataset Preparation CRN2 (2 Source)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

sensor_placement_matrix = zeros(M, N);
for i = 1:M
    sensor_placement_matrix(i, sensor_locations(i)+1) = 1;
end

numOfData = 3100000;
features = zeros(M*M, numOfData);
labels = zeros(2*N-1, numOfData);

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

    A = DOA.Array_Manifold(0.5, 0:N-1, doa);
    s = DOA.Source_Generate(K, L);
    Rs = (1 / L) * (s * s');
    R = A * Rs * A';

    labels(1:N, idx) = real(R(:, 1));
    labels(N+1:end, idx) = imag(R(2:end, 1));
    
    SNR_dB = 0 + rand * 30;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    A_ohm = sensor_placement_matrix * A;
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');

    re_R = real(R_ohm);
    im_R = imag(R_ohm);

    features(1:M, idx) = diag(re_R);
    ind = M + 1;
    for i = 2:M
        for j = 1:i-1
            features(ind:ind+1, idx) = [re_R(i, j); im_R(i, j)];
            ind = ind + 2;
        end
    end
end