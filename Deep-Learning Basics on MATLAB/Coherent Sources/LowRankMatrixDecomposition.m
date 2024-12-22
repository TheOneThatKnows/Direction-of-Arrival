%% Article Name: 
% DOA Estimation of Coherent Sources via Low-Rank Matrix Decomposition 

%% Notes
% Deep learning model can be used in order to reconstruct R

%% Initialize

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Define The Sensor Array

M = 4; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations_ca = sort([sensor_locations_1 sensor_locations_2(2:end)]);
sensor_locations_oca = [0 max(M,N)+sensor_locations_ca(2:end)];

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
% alpha = eye(K_coherent);
B = DOA.kronecker(conj(A), A);
P = alpha * alpha';
p = P(:);

Bp = B * p;
PHI = Rearrange_According_to_Sensor_Locations(DOA, Bp, sensor_locations_oca);

R = zeros(N);
for i = 1:N
    start_idx = N + 1 - i;
    end_idx = start_idx + N - 1;
    R(:, i) = PHI(start_idx:end_idx);
end

%%

B22 = Rearrange_According_to_Sensor_Locations(DOA, B, sensor_locations_oca);

%%

R_temp = reshape(Bp, M, M);
R2 = toeplitz(Virtual_Covariance_Column(DOA, R_temp, sensor_locations_oca)');

A2 = DOA.Array_Manifold(0:N-1, doa);
% A2([15 17 19 20], :) = 0;
% R_out = R_Toeplitz(A2 * P * A2', "full");
% R_out = A2 * P * A2';
R_out = A2 * diag(diag(P)) * A2';

%%

R_out_2 = A2 * P * A2';
R_out_2 = R_Toeplitz(R_out_2, "full");
% R_out_2 = R_out_2 / R_out_2(1, 1);
% R = R / R(1,1);
% 
% norm(R-R_out_2, "fro")

%%

uDOF = DOA.Uniform_Degrees_Of_Freedom(sensor_locations_oca);
N2 = (uDOF+1)*0.5;
A3 = DOA.Array_Manifold(0:N2-1, doa);

Rs_3 = ((A3' * A3) \ A3') * R(1:N2, 1:N2) * (A3 / (A3' * A3));
R_out_3 = A2 * diag(diag(Rs_3)) * A2';

%%

v3 = R_out_3(:, 1);
v3([15 17 19 20]) = 0;
R_out_3 = toeplitz(v3');

%%
% R_out_3 = R_out_3 / R_out_3(1, 1);
% R = R / R(1,1);

norm(R / R(1,1) - R_out_3 / R_out_3(1, 1), "fro")

%%

virtual_sensor_locations = 0:N-1;
virtual_sensor_locations([15 17 19 20]) = 0;
virtual_sensor_locations = sort(virtual_sensor_locations, "ascend");
virtual_sensor_locations = virtual_sensor_locations(5:end);

A4 = DOA.Array_Manifold(virtual_sensor_locations, doa);
v4 = R(:, 1);
v4 = v4(virtual_sensor_locations + 1);
R_in = toeplitz(v4');


Rs_4 = ((A4' * A4) \ A4') * R_in * (A4 / (A4' * A4));
R_out_4 = A2 * diag(diag(Rs_4)) * A2';

%%

figure; hold on
ss = DOA.MUSIC(K, R(1:N2, 1:N2), 0:N2-1, angle_spec);
plot(angle_spec, 10*log10(ss))
ss = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss))
ss = DOA.MUSIC(K, R_out_2, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss))
ss = DOA.MUSIC(K, R_out_3, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss))
ss = DOA.MUSIC(K, R_out_4, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss))
legend('R', 'R_out', 'R_out_2', 'R_out_3', 'R_out_4')
% legend('R_out', 'R_out_3')

%%

PHI3 = Rearrange_According_to_Sensor_Locations2(DOA, Bp, sensor_locations_oca);

R3 = zeros(N);
for i = 1:N
    start_idx = N + 1 - i;
    end_idx = start_idx + N - 1;
    R3(:, i) = PHI3(start_idx:end_idx);
end

%% Functions

% Sort and Discard Repeating Rows (Take Average) According to Sensor Locations In The Difference Coarray
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

% Sort and Discard Repeating Rows (Ignore) According to Sensor Locations In The Difference Coarray
function [X2, M_v] = Rearrange_According_to_Sensor_Locations2(DOA, X1, sensor_locations)
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
        if diff_coarray(idx2) == 0 || X2(idx2) ~= 0
            continue
        end

        idx1 = (i-1) * M + j;
        X2(idx2, :) = X1(idx1, :);
    end
end
end

% function R_out = SetLag2Zero(R, lags)
% M = size(R, 1);
% for i = 1:M
%     for j = 1:length(lags)
%         lag = lags(j);
%         R(i, i+lags) = 0;
%     end
% end
% end