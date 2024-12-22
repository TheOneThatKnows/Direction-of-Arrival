%% Define The Sensor Array

clear; clc; close all

M = 4; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations_ca = sort([sensor_locations_1 sensor_locations_2(2:end)]);
sensor_locations_oca = [0 max(M,N)+sensor_locations_ca(2:end)]

M1 = floor((M+N-1)/2);
M2 = ceil((M+N-1)/2);

ula1 = 0:M1-1;
ula2 = 0:M1+1:(M2-1)*(M1+1);
sensor_locations_na = [ula1 ula2+M1];
sensor_locations_na(M1+2:end) = sensor_locations_na(M1+2:end) + 1

%% Delete

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors
sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;
K_coherent = 2;
L = 70;

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A = DOA.Array_Manifold(sensor_locations, doa);
vars = ones(K-K_coherent+sign(K_coherent), 1);
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);

y = A * s + n;
Ry = (1 / L) * (y * y');

R_top = toeplitz(Virtual_Covariance_Column(DOA, Ry, sensor_locations));

%%

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;
K_coherent = 2;
L = 70;

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A = DOA.Array_Manifold(sensor_locations, doa);
vars = ones(K-K_coherent+sign(K_coherent), 1);
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);

y = A * s + n;
Ry = (1 / L) * (y * y');

R_top = R_Toeplitz(Ry, "full");

R_corr = zeros(M);
for i = 1:M
    for j = 1:M
        R_corr(i, j) = (1 / (2 * L - 1)) * (xcorr(y(i, :), y(j, :)) * xcorr(y(i, :), y(i, :))');
    end
end
R_ctop = R_Toeplitz(R_corr, "full");

ss1 = DOA.MUSIC(K, R_top, sensor_locations, angle_spec);
ss2 = DOA.MUSIC(K, R_ctop, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(ss1))
plot(angle_spec, 10*log10(ss2))
legend('TOP', 'CTOP')

% Ry_incoherent = A * diag(diag(Rs)) * A' + (1 / SNR) * eye(M)

% [eig_vecs, eig_vals] = eig(Ry_incoherent, "vector");
% [eig_vals, sorted_inds] = sort(eig_vals, "descend");
% eig_vecs = eig_vecs(:, sorted_inds)
% 
% eig_vecs(:, 1:K)' * eig_vecs(:, 1:K)

%% 

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:5; % ULA with 6 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A = DOA.Array_Manifold(sensor_locations, doa);

Rs = [1 0.8 0; 0.8 0.64 0; 0 0 1];

SNR_dB = 10;
SNR = 10^(SNR_dB / 10);
Ry = A * Rs * A' + (1 / SNR) * eye(M)
Ry_toeplitz = R_Toeplitz(Ry, "full")
Ry_incoherent = A * diag(diag(Rs)) * A' + (1 / SNR) * eye(M)

disp("Ry, Ry_toeplitz: ")
norm(Ry - Ry_toeplitz, "fro")

disp("Ry, Ry_incoherent: ")
norm(Ry - Ry_incoherent, "fro")

disp("Ry_toeplitz, Ry_incoherent: ")
norm(Ry_toeplitz - Ry_incoherent, "fro")

[eig_vecs, eig_vals] = eig(Ry, "vector")
[eig_vecs, eig_vals] = eig(Ry_toeplitz, "vector")
[eig_vecs, eig_vals] = eig(Ry_incoherent, "vector")

%% Initialization (ULA)

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:5; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

% doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
doa = [70 85 130];

vars = ones(K, 1);
vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

Ry = (1 / L) * (y * y');
R_fb = FBSS(Ry, 4);
R_toeplitz = R_Toeplitz(Ry, "full")

ss1 = DOA.MUSIC(K, Ry, 0:M-1, angle_spec);
ss2 = DOA.MUSIC(K, R_fb, 0:3, angle_spec);
ss3 = DOA.MUSIC(K, R_toeplitz, 0:M-1, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(ss1));
plot(angle_spec, 10*log10(ss2));
plot(angle_spec, 10*log10(ss3));
legend('R_y', 'R_{fb}', 'R_{toeplitz}');
%% 
Rs1 = (1 / L) * (s * s')
Rn1 = (1 / L) * (n * n');
Ry1 = (1 / L) * (y * y');

Rs2 = diag(diag(Rs1))
SNR = 10^(SNR_dB / 10);
Rn2 = (1 / SNR) * eye(M);
% Ry1_toeplitz = R_Toeplitz(Ry1, "full");
Ry1_2 = A * Rs1 * A' + Rn1;
Ry1_toeplitz = R_Toeplitz(Ry1_2, "full");
Ry2 = A * Rs2 * A' + Rn1;

disp("Ry1_2, Ry1_toeplitz: ")
norm(Ry1_2 - Ry1_toeplitz, "fro")

disp("Ry2, Ry1_toeplitz: ")
norm(Ry2 - Ry1_toeplitz, "fro")

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

% s = DOA.Source_Generate_old(K, L);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 20;
% n = DOA.Noise_Generate(SNR_dB, M, L);
% y = A * s + n;

Rs = [1 0 0; 0 1 1; 0 1 1];
SNR = 10^(SNR_dB / 10);
Ry = A * Rs * A' + (1 / SNR) * eye(M);

Rs_out = R_Toeplitz(Rs, "full");
Ry_out = R_Toeplitz(Ry, "full");
Ry_out_2 = A * Rs_out * A' + (1 / SNR) * eye(M);

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

K = 3;          % # of sources
K_coherent = 3;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

% s = DOA.Source_Generate_old(K, L);
% shuffledIndices = randperm(K);
% s = s(shuffledIndices, :);

vars = ones(K, 1);
vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
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

% A = DOA.Array_Manifold(0:N-1, doa);
% R = A * A';
% spatial_spectrum_3 = DOA.MUSIC(K, R, 0:N-1, angle_spec);

R = toeplitz(Virtual_Covariance_Column(DOA, R_ohm, sensor_locations)');
spatial_spectrum_3 = DOA.MUSIC(K, R, 0:N-1, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_4))
plot(angle_spec, 10*log10(spatial_spectrum_1))
plot(angle_spec, 10*log10(spatial_spectrum_2))
plot(angle_spec, 10*log10(spatial_spectrum_3))
legend('R_ohm', 'SS', 'Toeplitz', 'VCC')

% figure; hold on
% plot(angle_spec, 10*log10(spatial_spectrum_1))




%%

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
sensor_locations = 0:4; % MRA with 5 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 4;
K_coherent = 2;

alpha = 1 - (0:0.1:0.1*(K_coherent-1)).';
Rs = eye(K);
Rs(1:K_coherent, 1:K_coherent) = alpha * alpha';

sum1 = zeros(M, 1);
sum2 = zeros(M, 1);
EPOCHS = 5000;
for i = 1:EPOCHS
    doa = DOA.DOA_Generate(K, 30, 150, 1);
    A_ohm = DOA.Array_Manifold(sensor_locations, doa);
    R = A_ohm * Rs * A_ohm';
    rank(R)

    [~, eig_vals_1] = eig(R, "vector");

    % R_toeplitz = Sparse_Toeplitz(DOA, R, sensor_locations);
    R_toeplitz = R_Toeplitz(R, "half");
    rank(R_toeplitz)

    [~, eig_vals_2] = eig(R_toeplitz, "vector");

    sum1 = sum1 + sort(eig_vals_1, "descend");
    sum2 = sum2 + sort(eig_vals_2, "descend");
end

sum1 = sum1 / EPOCHS
sum2 = sum2 / EPOCHS

% [sort(eig_vals_1, "descend") sort(eig_vals_2, "descend")]

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [1 2 3 4 9 13] - 1; % NA_v2 with 6 sensors
sensor_locations = [0 1 4 7 9]; % NA_v2 with 6 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 4;          % # of sources
K_coherent = 3;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

vars = [1 1];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
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

R_sparse_toeplitz = Sparse_Toeplitz(DOA, R_ohm2_noiseless, sensor_locations);
disp("Rank(R_sparse_toeplitz): " + rank(R_sparse_toeplitz))
[~, eig_vals] =  eig(R_sparse_toeplitz, "vector");
eig_vals_R_sparse_toeplitz = sort(eig_vals, "descend");

R_ss2 = Spatial_Smoothing(DOA, sensor_locations, R_ohm2_noiseless);
disp("Rank(R_ss2): " + rank(R_ss2))
[~, eig_vals] =  eig(R_ss2, "vector");
eig_vals_R_ss2 = sort(eig_vals, "descend");

% disp([[eig_vals_ohm2; zeros(N-M, 1)] [eig_vals_ohm2_noiseless; zeros(N-M, 1)] eig_vals_R])
disp([[eig_vals_ohm2; zeros(N-M, 1)] eig_vals_R eig_vals_R_ss2])

%%

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [1 2 3 4 9 13] - 1; % NA_v2 with 6 sensors
sensor_locations = [0 1 4 7 9]; % NA_v2 with 6 sensors

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
% R_ss = Spatial_Smoothing(DOA, sensor_locations, R_ohm);
R_sparse_toeplitz = Sparse_Toeplitz(DOA, R_ohm, sensor_locations);

spatial_spectrum_1 = DOA.MUSIC(K, R_ohm, sensor_locations, angle_spec);
spatial_spectrum_2 = DOA.SS_MUSIC(K, R_ohm, sensor_locations, angle_spec);
spatial_spectrum_3 = DOA.MUSIC(K, R_sparse_toeplitz, sensor_locations, angle_spec);
spatial_spectrum_4 = DOA.SS_MUSIC(K, R_sparse_toeplitz, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_1))
plot(angle_spec, 10*log10(spatial_spectrum_2))
plot(angle_spec, 10*log10(spatial_spectrum_3))
plot(angle_spec, 10*log10(spatial_spectrum_4))
legend('MUSIC', 'SS-MUSIC', 'Sparse-Toeplitz', 'SS-MUSIC-2')

%%

spatial_spectrum_4 = DOA.DML(R_ohm, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10*log10(spatial_spectrum_4))

