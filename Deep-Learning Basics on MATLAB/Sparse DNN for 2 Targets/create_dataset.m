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
features = zeros(N*N, numOfData);

K = 2;          % # of sources
L = 100;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
labels = zeros((phi_max - phi_min)/delta_phi + 1, numOfData);
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
    SNR_dB = 30 * rand;     % min = 0 dB, max = 30 dB
    n = DOA.Noise_Generate(SNR_dB, M, L);
    y = A * s + n;
    R = (1 / L) * (y * y');

    z = R(:);
    z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
    R_z1 = zeros(N);
    for i = 1:N
        z1_i = z1(i:i + N - 1);
        R_z1 = R_z1 + (1 / N) * (z1_i * z1_i');
    end

    re_R_z1 = real(R_z1);
    im_R_z1 = imag(R_z1);

    features(1:N, idx) = diag(re_R_z1);
    ind = N + 1;
    for i = 2:N
        for j = 1:i-1
            features(ind:ind+1, idx) = [re_R_z1(i, j); im_R_z1(i, j)];
            ind = ind + 2;
        end
    end
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

function h5_Format(var, filename, dataset)
filename = filename + ".h5";
h5create(filename, dataset, size(var));
h5write(filename, dataset, var);
h5disp(filename);
end