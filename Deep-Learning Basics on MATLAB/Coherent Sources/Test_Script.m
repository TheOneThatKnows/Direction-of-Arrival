%% Initialization (ULA)

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

s = DOA.Source_Generate_old(K, L);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 20;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

Ry = (1 / L) * (y * y');

spatial_spectrum_1 = DOA.MUSIC(K, Ry, sensor_locations, angle_spec);

R_out = R_Toeplitz(Ry, "half");
spatial_spectrum_2 = DOA.MUSIC(K, R_out, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_1))
plot(angle_spec, 10*log10(spatial_spectrum_2))
legend('Ry', 'Toeplitz')

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_1))









%% Initialization (Old Source Generate)

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [1 2 3 4 9 13] - 1; % NA_v2 with 6 sensors
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 2;          % # of sources
K_coherent = 0;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

s = DOA.Source_Generate_old(K, L);
% shuffledIndices = randperm(K);
% s = s(shuffledIndices, :);
A_ohm = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 20;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;

R_ohm = (1 / L) * (y * y');

spatial_spectrum_4 = DOA.MUSIC(K, R_ohm, sensor_locations, angle_spec);

% this one is just to calculate the rank
Rs = round((1 / L) * abs(s * s'));
SNR = 10^(SNR_dB / 10);
R_ohm2_noiseless = A_ohm * Rs * A_ohm';
R_ohm2 = R_ohm2_noiseless + (1 / SNR) * eye(M);

R_ss = Spatial_Smoothing(DOA, sensor_locations, R_ohm);
% spatial_spectrum_1 = DOA.MUSIC(K, R_ss, 0:N-1, angle_spec);
[spatial_spectrum_1, R_ss2] = DOA.SS_MUSIC(K, R_ohm, sensor_locations, angle_spec);

R_out = R_Toeplitz(R_ss, "half");
spatial_spectrum_2 = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);

A = DOA.Array_Manifold(0:N-1, doa);
R = A * A';
spatial_spectrum_3 = DOA.MUSIC(K, R, 0:N-1, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_4))
plot(angle_spec, 10*log10(spatial_spectrum_1))
plot(angle_spec, 10*log10(spatial_spectrum_2))
plot(angle_spec, 10*log10(spatial_spectrum_3))
legend('R_ohm', 'SS', 'Toeplitz', 'ULA')

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_1))








%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = [1 2 3 4 9 13] - 1; % NA_v2 with 6 sensors

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
s_old = s;
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
A_ohm = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;

R_ohm = (1 / L) * (y * y');

% this one is just to calculate the rank
Rs = round((1 / L) * abs(s * s'));
SNR = 10^(SNR_dB / 10);
R_ohm2_noiseless = A_ohm * Rs * A_ohm';
R_ohm2 = R_ohm2_noiseless + (1 / SNR) * eye(M);

R_ss = Spatial_Smoothing(DOA, sensor_locations, R_ohm);
% spatial_spectrum_1 = DOA.MUSIC(K, R_ss, 0:N-1, angle_spec);
spatial_spectrum_1 = DOA.SS_MUSIC(K, R_ohm, sensor_locations, angle_spec);

R_out = R_Toeplitz(R_ss, "half");
spatial_spectrum_2 = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);

A = DOA.Array_Manifold(0:N-1, doa);
R = A * A';
spatial_spectrum_3 = DOA.MUSIC(K, R, 0:N-1, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_1))
plot(angle_spec, 10*log10(spatial_spectrum_2))
plot(angle_spec, 10*log10(spatial_spectrum_3))
legend('SS', 'Toeplitz', 'ULA')

%%

disp(Rs)
disp("Rank(Rs): " + rank(Rs))

disp("Rank(R_ohm2): " + rank(R_ohm2))
[~, eig_vals] =  eig(R_ohm2, "vector");
eig_vals_ohm2 = sort(eig_vals, "descend");

disp("Rank(R_ohm2_noiseless): " + rank(R_ohm2_noiseless))
[~, eig_vals] =  eig(R_ohm2_noiseless, "vector");
eig_vals_ohm2_noiseless = sort(eig_vals, "descend");

disp("Rank(R): " + rank(R))
[~, eig_vals] =  eig(R, "vector");
eig_vals_R = sort(eig_vals, "descend");

disp([[eig_vals_ohm2; zeros(N-M, 1)] [eig_vals_ohm2_noiseless; zeros(N-M, 1)] eig_vals_R])

%%

spatial_spectrum_4 = DOA.DML(R_ohm, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_4))