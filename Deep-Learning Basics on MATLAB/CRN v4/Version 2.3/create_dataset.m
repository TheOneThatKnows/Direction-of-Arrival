%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

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
% numOfData = 200;

delta_phi = 3;
phi_min = 30;
phi_max = 150;
Q = (phi_max - phi_min) / delta_phi + 1;

features = zeros(M, M, 3, numOfData);
labels = zeros(Q, numOfData);

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

    for i = 1:K
        c = floor((doa(i) - phi_min) / delta_phi);

        phi_a = phi_min + c * delta_phi;
        phi_b = phi_a + delta_phi;

        percentages = [phi_a phi_b; 1 1] \ [doa(i)*100; 100];

        idx1 = c + 1;
        idx2 = c + 2;
        labels(idx1:idx2, idx) = labels(idx1:idx2, idx) + percentages;
    end

    A_ohm = DOA.Array_Manifold(0.5, sensor_locations, doa);
    L = 70; % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A_ohm * s + n;
    R_ohm = (1 / L) * (y * y');
    normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

    features(:, :, 1, idx) = real(normalized_R_ohm);
    features(:, :, 2, idx) = imag(normalized_R_ohm);
    features(:, :, 3, idx) = angle(R_ohm) / pi;
end