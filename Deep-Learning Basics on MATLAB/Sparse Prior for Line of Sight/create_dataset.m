%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();

%% Sensor Properties

% sensor_locations = [0 1 2 6 10 13];

%% Dataset Preparation I

sensor_locations = 0:7;
M = length(sensor_locations);
M2 = M * M;

numOfData = 6000;
features = zeros((M-1)*M, numOfData);
labels = zeros(M2*M2, numOfData);

L = 200;
K = 2;
delta_phi = 1;
phi_min = 80;
phi_max = 100;
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

    % original case
    sensor_locations = 0:7;
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    s = DOA.Source_Generate(K, L);
    SNR_dB = -10 + rand * 30;
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
    r = (r - mean(r)) / std(r);

    features(:, idx) = r;

    % desired case
    sensor_locations = 0:M2-1;
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    y = A * s;
    R = (1 / L) * (y * y');

    re_R = real(R);
    im_R = imag(R);
    re_r = zeros((M2 + 1) * M2 / 2, 1);
    im_r = zeros((M2 - 1) * M2 / 2, 1);

    idx2 = 0;
    for i = 1:M2
        idx1 = idx2 + 1;
        idx2 = idx1 + M2 - i;
        re_r(idx1:idx2) = re_R(i, i:M2);
    end

    idx2 = 0;
    for i = 1:M2-1
        idx1 = idx2 + 1;
        idx2 = idx1 + M2 - i - 1;
        im_r(idx1:idx2) = im_R(i, i+1:M2);
    end
    r = [re_r; im_r];
    r = (r - mean(r)) / std(r);

    labels(:, idx) = r;
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