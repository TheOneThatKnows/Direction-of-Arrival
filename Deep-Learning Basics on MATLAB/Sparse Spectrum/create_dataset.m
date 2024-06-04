%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors
% sensor_locations = [0 1 2 6 10 13]; % SLA with 6 sensors

%% Dataset Preparation

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = 1;
phi_min = 30;
phi_max = 150;
angles = phi_min:delta_phi:phi_max;
Q = length(angles);

numOfData = 43200;
numOfData = 200;
features = zeros(Q, 1, numOfData);
labels = zeros(Q, numOfData);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angles);
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
    rounded_doa = round(doa);
    inds = (rounded_doa - phi_min) / delta_phi + 1;
    labels(inds, idx) = 1;

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    % L = round(randn * 30 + 500); % # of snapshots
    L = 70;
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');
    
    z = R(:);
    A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);

    feature = abs(A1' * z);
    features(:, 1, idx) = (feature - mean(feature)) / var(feature);
end

%% 

figure; hold on;
plot(angles, features(:, 1, numOfData));

[z1, M_v] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
A2 = DOA.Array_Manifold(0.5, -M_v+1:M_v-1, angles);

feature = abs(A2' * z1);
feature = (feature - mean(feature)) / var(feature);

plot(angles, feature);

legend('z', 'z1')

%% 

pred = predict(net, feature);

figure; hold on;
plot(angles, feature);
plot(angles, pred);
xlabel('angles');
legend('feature', 'output');
grid on;

msg = "Actual DOAs:";
for i = 1:K
    msg = msg + " " + doa(i);
end
title('Feature Spectrum and Output Spectrum (' + msg + ')');

%% Dataset Preparation

% ilk olarak 30 ile 150 arası feature ı kullan
% sonra 87 ile 90 arası feature kullan
% birleştir

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = [1 0.025];
phi_min = [30 87];
phi_max = [150 90];
angles = phi_min(2):delta_phi(2):phi_max(2);
Q = length(angles);

numOfData = 43200;
features = zeros(Q, 1, numOfData);
labels = zeros(Q, numOfData);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angles);
for idx = 1:numOfData
    doa = zeros(1, K);
    doa(1) = phi_min(2) + rand * (phi_max(2) - phi_min(2));
    i = 2;
    while true
        temp_angle = phi_min(1) + rand * (phi_max(1) - phi_min(1));
        if temp_angle > phi_min(2) && temp_angle < phi_max(2)
            i = i - 1;
        else
            doa(i) = temp_angle;
        end
        if i == K
            break
        end
        i = i + 1;
    end
    rounded_doa = doa;
    for i = 1:K
        remainder = rem(rounded_doa(i), delta_phi(2));
        if remainder < delta_phi(2) / 2
            rounded_doa(i) = rounded_doa(i) - remainder;
        else
            rounded_doa(i) = rounded_doa(i) + delta_phi(2) - remainder;
        end
    end
    ind = round((rounded_doa(1) - phi_min(2)) / delta_phi(2) + 1);
    labels(ind, idx) = 1;

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    % L = round(randn * 30 + 500); % # of snapshots
    L = 70;
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');
    
    z = R(:);
    A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);

    feature = abs(A1' * z);
    features(:, 1, idx) = (feature - mean(feature)) / var(feature);
end

%%

figure; hold on;
plot(angles, features(:, :, end));
plot(angles, labels(:, end));

%% Dataset Preparation

% ilk olarak 30 ile 150 arası feature ı kullan
% sonra 87 ile 90 arası feature kullan
% birleştir

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = [0.25 0.25];
phi_min = [30 87];
phi_max = [150 90];
angles = phi_min(1):delta_phi(2):phi_max(1);
Q = length(angles);
Q2 = (phi_max(2) - phi_min(2)) / delta_phi(2) + 1;

numOfData = 43200;
features = zeros(Q, 1, numOfData);
labels = zeros(Q2, numOfData);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angles);
for idx = 1:numOfData
    doa = zeros(1, K);
    doa(1) = phi_min(2) + rand * (phi_max(2) - phi_min(2));
    i = 2;
    while true
        temp_angle = phi_min(1) + rand * (phi_max(1) - phi_min(1));
        if temp_angle > phi_min(2) && temp_angle < phi_max(2)
            i = i - 1;
        else
            doa(i) = temp_angle;
        end
        if i == K
            break
        end
        i = i + 1;
    end
    rounded_doa = doa;
    for i = 1:K
        remainder = rem(rounded_doa(i), delta_phi(2));
        if remainder < delta_phi(2) / 2
            rounded_doa(i) = rounded_doa(i) - remainder;
        else
            rounded_doa(i) = rounded_doa(i) + delta_phi(2) - remainder;
        end
    end
    ind = round((rounded_doa(1) - phi_min(2)) / delta_phi(2) + 1);
    labels(ind, idx) = 1;

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    % L = round(randn * 30 + 500); % # of snapshots
    L = 70;
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');
    
    z = R(:);
    A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);

    feature = abs(A1' * z);
    features(:, 1, idx) = (feature - mean(feature)) / var(feature);
end

%% Dataset Preparation

% ilk olarak 30 ile 150 arası feature ı kullan
% sonra 87 ile 90 arası feature kullan
% birleştir

M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;  % # of sources

delta_phi = [0.25 0.25];
phi_min = [30 87];
phi_max = [150 90];
angles = phi_min(1):delta_phi(2):phi_max(1);
Q = length(angles);

numOfData = 43200;
features = zeros(Q, 1, numOfData);
labels = zeros(Q, numOfData);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angles);
for idx = 1:numOfData
    doa = zeros(1, K);
    doa(1) = phi_min(2) + rand * (phi_max(2) - phi_min(2));
    i = 2;
    while true
        temp_angle = phi_min(1) + rand * (phi_max(1) - phi_min(1));
        if temp_angle > phi_min(2) && temp_angle < phi_max(2)
            i = i - 1;
        else
            doa(i) = temp_angle;
        end
        if i == K
            break
        end
        i = i + 1;
    end
    rounded_doa = doa;
    for i = 1:K
        remainder = rem(rounded_doa(i), delta_phi(2));
        if remainder < delta_phi(2) / 2
            rounded_doa(i) = rounded_doa(i) - remainder;
        else
            rounded_doa(i) = rounded_doa(i) + delta_phi(2) - remainder;
        end
    end
    ind = round((rounded_doa(1) - phi_min(1)) / delta_phi(1) + 1);
    labels(ind, idx) = 1;

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    % L = round(randn * 30 + 500); % # of snapshots
    L = 70;
    s = DOA.Source_Generate(K, L);
    SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');
    
    z = R(:);
    A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);

    feature = abs(A1' * z);
    features(:, 1, idx) = (feature - mean(feature)) / var(feature);
end