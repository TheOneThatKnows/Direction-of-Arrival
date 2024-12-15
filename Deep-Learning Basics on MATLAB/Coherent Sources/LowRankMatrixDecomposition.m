%% Article Name: 
% DOA Estimation of Coherent Sources via Low-Rank Matrix Decomposition 

%% Initialize

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Define The Sensor Array

M = 4; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations_ca = sort([sensor_locations_1 sensor_locations_2(2:end)]);
sensor_locations_oca = [0 max(M,N)+sensor_locations_ca(1:end)];

M = length(sensor_locations_oca);
N = sensor_locations_oca(M) + 1;

%% Simulate The Environment

K = 3;          % # of sources
K_coherent = K;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

vars = ones(K - K_coherent + 1, 1);
% vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);
A = DOA.Array_Manifold(sensor_locations_oca, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

%% Start The Algorithm 

Ry = (1 / L) * (y * y');

r = Ry(:);
z = Rearrange_According_to_Sensor_Locations(DOA, r, sensor_locations_oca);

R = zeros(N);
for i = 1:N
    start_idx = N + 1 - i;
    end_idx = start_idx + N - 1;
    R(:, i) = z(start_idx:end_idx);
end

% The optimization starts here

%%

alpha = [1; (0.75+0.25*rand(K_coherent-1, 1)) .* exp(1i * pi * rand(K_coherent-1, 1))];
B = DOA.kronecker(conj(A), A);
p = alpha * alpha';
p = p(:);

PHI = Rearrange_According_to_Sensor_Locations(DOA, B * p, sensor_locations_oca);

R = zeros(N);
for i = 1:N
    start_idx = N + 1 - i;
    end_idx = start_idx + N - 1;
    R(:, i) = z(start_idx:end_idx);
end

%% Functions

% Sort and Discard Repeating Rows According to Sensor Locations In The Difference Coarray
function [X2, M_v] = Rearrange_According_to_Sensor_Locations(DOA, X1, sensor_locations)
diff_coarray = DOA.Diff_Coarray(sensor_locations);
noOfRows = length(diff_coarray);
M_v = 0.5 * (noOfRows + 1);
M = length(sensor_locations);

[~, noOfCols] = size(X1);
X2 = zeros(noOfRows, noOfCols);
for i = 1:length(sensor_locations)
    for j = 1:length(sensor_locations)
        diff = -sensor_locations(i) + sensor_locations(j);
        idx2 = M_v + diff;
        if diff_coarray(idx2) == 0
            continue
        end

        idx1 = (i-1) * M + j;
        X2(idx2, :) = X2(idx2, :) + (1/diff_coarray(idx2)) * X1(idx1, :);
    end
end
end