%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Sensor Properties

% sensor_locations = [0 1 2 6 10 13];
sensor_locations = 0:63;

%% Dataset Preparation I (This one is not used)

numOfData = sensor_locations(M) * 10000;
features = zeros((M-1)*M, numOfData);

delta_phi = 1;
phi_max = 180;
phi_min = 0;
labels = zeros((phi_max - phi_min)/delta_phi + 1, numOfData);
for K = 1:sensor_locations(M)
    for idx = 1:10000
        doa = (-2 * delta_phi) * ones(1, K);
        i = 1;
        while true
            temp_angle = rand * 180;
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
            idx1 = floor(doa(i) / delta_phi) + 1;
            idx2 = ceil(doa(i) / delta_phi) + 1;
            x0 = (idx1 - 1) * delta_phi;
            x1 = (idx2 - 1) * delta_phi;
            z = [1 1; x0 x1] \ [1; doa(i)];
            labels(idx1:idx2, (K-1)*10000 + idx) = z;
        end

        A = DOA.Array_Manifold(coef, sensor_locations, doa);
        s = DOA.Source_Generate(K, snapshots);
        n = DOA.Noise_Generate(SNR_dB, M, snapshots);
        y = A * s + n;
        R = (1 / snapshots) * (y * y');

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

        features(:, (K-1)*10000 + idx) = r;
    end
end

%% Dataset Preparation II
% A Novel Tree Model-based DNN to Achieve a High-Resolution DOA Estimation 
% via Massive MIMO receive array

numOfData = 100000;
[~, features_1, labels_1] = TreeDataset(DOA, numOfData, sensor_locations, [30 149 10]);
[~, features_2_1, labels_2_1] = TreeDataset(DOA, numOfData, sensor_locations, [30 39 1]);
[~, features_2_2, labels_2_2] = TreeDataset(DOA, numOfData, sensor_locations, [40 49 1]);
[~, features_2_3, labels_2_3] = TreeDataset(DOA, numOfData, sensor_locations, [50 59 1]);
[~, features_2_4, labels_2_4] = TreeDataset(DOA, numOfData, sensor_locations, [60 69 1]);
[~, features_2_5, labels_2_5] = TreeDataset(DOA, numOfData, sensor_locations, [70 79 1]);
[~, features_2_6, labels_2_6] = TreeDataset(DOA, numOfData, sensor_locations, [80 89 1]);
[~, features_2_7, labels_2_7] = TreeDataset(DOA, numOfData, sensor_locations, [90 99 1]);
[~, features_2_8, labels_2_8] = TreeDataset(DOA, numOfData, sensor_locations, [100 109 1]);
[~, features_2_9, labels_2_9] = TreeDataset(DOA, numOfData, sensor_locations, [110 119 1]);
[~, features_2_10, labels_2_10] = TreeDataset(DOA, numOfData, sensor_locations, [120 129 1]);
[~, features_2_11, labels_2_11] = TreeDataset(DOA, numOfData, sensor_locations, [130 139 1]);
[~, features_2_12, labels_2_12] = TreeDataset(DOA, numOfData, sensor_locations, [140 149 1]);

%% h5 FORMATTING

filename = "features_1.h5";
h5create(filename,"/DS1",[4032 100000]);
h5write(filename,"/DS1", features);
h5disp(filename);

filename = "labels_1.h5";
h5create(filename,"/DS1",[12 100000]);
h5write(filename,"/DS1", labels_1);
h5disp(filename);

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

function [angles, features, labels] = TreeDataset(DOA, numOfData, sensor_locations, phi_specifics)
angles = zeros(1, numOfData);

M = length(sensor_locations);
features = zeros((M-1)*M, numOfData);

K = 1; % # of sources is 1
phi_min = phi_specifics(1);
phi_max = phi_specifics(2);
delta_phi = phi_specifics(3);
labels = zeros((phi_max - phi_min + 1)/delta_phi, numOfData);

for idx = 1:numOfData
    doa = rand * (phi_max - phi_min) + phi_min;
    angles(idx) = doa;
    rounded_doa = round(doa);

    l_h = ceil((rounded_doa - phi_min + 1) / delta_phi);
    labels(l_h, idx) = 1;

    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    SNR_dB = randn + 10;
    SNR = 10^(SNR_dB / 10);
    
    R = A * A' + (1 / SNR) * eye(M);

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
end
end