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

numOfData = 1100000;
features = zeros(2*N-1, numOfData);

L = 100;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
labels = zeros(1, numOfData);
for idx = 1:numOfData
    K = randi(N-1);     % # of sources
    labels(idx) = K;
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

    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = 30 * rand;     % min = 0 dB, max = 30 dB
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    z = R(:);
    z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
    z2 = z1(1:N-1);
    features(:, idx) = [real(z2); z1(N); imag(z2)];
end

%% Dataset Preparation II

numOfData = 1100000;
features = zeros(M*M, numOfData);

L = 100;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
labels = zeros(1, numOfData);
for idx = 1:numOfData
    K = randi(N-1);     % # of sources
    labels(idx) = K;
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

    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = 30 * rand;     % min = 0 dB, max = 30 dB
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    re_R = real(R);
    im_R = imag(R);

    r(1:M) = diag(re_R);
    ind = M + 1;
    for i = 2:M
        for j = 1:i-1
            r(ind:ind+1) = [re_R(i, j); im_R(i, j)];
            ind = ind + 2;
        end
    end
    features(:, idx) = r;
end