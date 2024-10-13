%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

%% Dataset Preparation CRN2 (2 Source)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = 1200000;

features = zeros(M, M, 3, numOfData);
labels = zeros(2 * N - 1, numOfData);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
for idx = 1:numOfData
    K = 2;
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

    A = DOA.Array_Manifold(0.5, 0:N-1, doa);
    L = 70; % # of snapshots
    s = DOA.Source_Generate(K, L);
    Rs = (1 / L) * (s * s');
    R = A * Rs * A';

    % converting R into a toeplitz matrix
    R_T = zeros(N);
    for m = -N+1:N-1
        J = zeros(N);
        J_temp = eye(N - abs(m));
        J(1:N-abs(m), abs(m)+1:N) = J_temp;
        if m < 0
            J = J.';
        end
        R_T = R_T + 1 / (N - abs(m)) * trace(R * J) * J.';
    end

    label = zeros(2 * N - 1, 1);
    label(1:N) = real(R_T(:, 1));
    imag_R_T = imag(R_T);
    label(N+1:end) = imag_R_T(2:N, 1);
    labels(:, idx) = label / label(1);

    A_ohm = DOA.Array_Manifold(0.5, sensor_locations, doa);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');
    normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

    features(:, :, 1, idx) = real(normalized_R_ohm);
    features(:, :, 2, idx) = imag(normalized_R_ohm);
    features(:, :, 3, idx) = angle(R_ohm) / pi;
end