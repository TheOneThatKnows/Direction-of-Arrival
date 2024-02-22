%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();

%% Sensor Properties

sensor_locations = [1 2 5 7] - 1; % SLA with 4 sensors
% sensor_locations = [0 1 2 6 10 13]; % SLA with 6 sensors

%% Dataset Preparation TDN

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = 600000;
features = zeros(M, M, 2, numOfData);
labels = zeros(1, numOfData);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
for idx = 1:numOfData
    K = randi(N-1);
    labels(idx) = K;
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

    % original case
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    L = round(randn * 30 + 500); % # of snapshots
    s = DOA.Source_Generate(K, L);
    SNR_dB = -15 + rand * 30;
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    features(:, :, 1, idx) = real(R);
    features(:, :, 2, idx) = imag(R);
end

%% h5 FORMATTING

h5_Format(features_1, "features_1", "/DS1");
h5_Format(labels_1, "labels_1", "/DS1");

%% JSON FORMATTING (did not work)

% Encode the data
encodedFeatures = jsonencode(features);
encodedLabels1 = jsonencode(labels_1);
encodedLabels2 = jsonencode(labels_2);

% Write to a JSON file
filename = 'features.json';
fid = fopen(filename, 'w');
fprintf(fid, '%s', encodedFeatures);
fclose(fid);

disp(['Data saved to ' filename]);

% Repeat the same process for labels_1
filename = 'labels1.json';
fid = fopen(filename, 'w');
fprintf(fid, '%s', encodedLabels1);
fclose(fid);

disp(['Data saved to ' filename]);

%% CSV FORMATTING (did not work)

filename = 'features';
writematrix(features, filename)
disp(['Data saved to ' filename]);

filename = 'labels1';
writematrix(labels_1, filename)
disp(['Data saved to ' filename]);

%% Functions

function r = Concatenate_Covariance_Vector(DOA, sensor_locations, doa, K, L, SNR_dB, M)
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
s = DOA.Source_Generate(K, L);
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;
R = (1 / L) * (y * y');

re_R = real(R);
im_R = imag(R);
re_r = zeros((M-1) * M / 2, 1);
im_r = zeros((M-1) * M / 2, 1);

idx2 = 0;
for i = 1:M-1
    idx1 = idx2 + 1;
    idx2 = idx1 + M - i - 1;
    re_r(idx1:idx2) = re_R(i, i+1:M);
    im_r(idx1:idx2) = im_R(i, i+1:M);
end
r = [re_r; im_r];
end

function h5_Format(var, filename, dataset)
filename = filename + ".h5";
h5create(filename, dataset, size(var));
h5write(filename, dataset, var);
h5disp(filename);
end