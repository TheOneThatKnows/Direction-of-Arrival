%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Sensor Properties

sensor_locations = [0 1 2 6 10 13];
M = length(sensor_locations);

% sensor_locations = DOA.Sensor_Locations(input)                  % sensor locations
sensor_placement = DOA.Sensor_Placement(sensor_locations)       % sensor places (in order to visualize the array)
diff_coarray = DOA.Diff_Coarray(sensor_placement)               % difference coarray
uDOF = DOA.Uniform_Degrees_Of_Freedom(diff_coarray)             % uniform degrees of freedom
LuDOF = DOA.One_Side_Uniform_Degrees_Of_Freedom(diff_coarray)   % one side uniform degrees of freedom

snapshots = 200;                                            % # of snapshots
SNR_dB = 10;                                                % signal to noise ratio in decibels
C = DOA.Mutual_Coupling(0, 0.1, M, sensor_locations);     % mutual coupling

%% Dataset Preparation

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