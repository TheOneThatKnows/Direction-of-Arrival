%% Article

% Title: 
% Direction-of-Arrival Estimation of Mixed Coherent and Uncorrelated Signals 

% Link: 
% https://ieeexplore.ieee.org/abstract/document/10640208?casa_token=8V0l9DteRzYAAAAA:QzEjv5QNOgmTkZ3-I4KgSoOrpftMW9hfe6YO8UcnNKHG59qPupgeLd59_G4e9rG0Kj5Ig-fw6w

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:6; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

s = [DOA.Source_Generate(K-K_coherent, L); 
    DOA.Coherent_Source_Generate(K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');

spec = DOA.MUSIC(K, R, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10 * log10(spec));

R_out = R_Toeplitz(R, "half");

spec = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
plot(angle_spec, 10 * log10(spec));
title_text = "DOAs: [";
for k = 1:K
    title_text = title_text + " " + doa(k);
end
title_text = title_text + " ]";
title(title_text)
legend('1', '2')

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
% sensor_locations = [1 2 3 4 8 12] - 1; % NA with 6 sensors
sensor_locations = [1 2 3 4 9 13] - 1; % NA v2 with 6 sensors
% sensor_locations = [1 2 3 4 5 11 16] - 1; % NA v2 with 8 sensors
% sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

s = [DOA.Source_Generate(K-K_coherent, L); 
    DOA.Coherent_Source_Generate(K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');

spec = DOA.MUSIC(K, R, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10 * log10(spec));

R = Spatial_Smoothing(DOA, sensor_locations, R);

R_out = R_Toeplitz(R, "half");

spec = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
plot(angle_spec, 10 * log10(spec));
title_text = "DOAs: [";
for k = 1:K
    title_text = title_text + " " + doa(k);
end
title_text = title_text + " ]";
title(title_text)
legend('1', '2')

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 10)

s = [DOA.Source_Generate(K-K_coherent, L); 
    DOA.Coherent_Source_Generate(K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
SNR_dB = 10;

% ULA
sensor_locations = 0:6; % ULA with 7 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

A = DOA.Array_Manifold(sensor_locations, doa);
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');
spec_ula_og = DOA.MUSIC(K, R, sensor_locations, angle_spec);

R_ula_full = R_Toeplitz(R, "full");
spec_ula_full = DOA.MUSIC(K, R_ula_full, sensor_locations, angle_spec);

R_ula_algo = R_Toeplitz(R, "half");
spec_ula_algo = DOA.MUSIC(K, R_ula_algo, sensor_locations, angle_spec);

% Sparse

sensor_locations = [1 2 3 4 9 13 17] - 1; % NA v2 with 7 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

A = DOA.Array_Manifold(sensor_locations, doa);
% n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');
spec_sparse_og = DOA.MUSIC(K, R, sensor_locations, angle_spec);

R = Spatial_Smoothing(DOA, sensor_locations, R);
spec_ss_sparse = DOA.MUSIC(K, R, 0:N-1, angle_spec);

R_sparse_algo = R_Toeplitz(R, "half");
spec_sparse_algo = DOA.MUSIC(K, R_sparse_algo, 0:N-1, angle_spec);

% Plot
figure; hold on
plot(angle_spec, 10 * log10(spec_ula_og));
plot(angle_spec, 10 * log10(spec_ula_full));
plot(angle_spec, 10 * log10(spec_ula_algo));
plot(angle_spec, 10 * log10(spec_sparse_og));
plot(angle_spec, 10 * log10(spec_ss_sparse));
plot(angle_spec, 10 * log10(spec_sparse_algo));
title_text = "DOAs: [";
for k = 1:K
    title_text = title_text + " " + doa(k);
end
title_text = title_text + " ]";
title(title_text)
legend('ULA Og', 'ULA full', 'ULA + Algo', 'Sparse Og', 'SS-Sparse', 'Sparse + Algo')