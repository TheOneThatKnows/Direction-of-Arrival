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

K = 2;          % # of sources
K_coherent = 2;
L = 200;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

% doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
doa = [80 90];
% doa = [50 79.2 90];

vars = ones(K - K_coherent + 1, 1);
% vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
alpha = [1; exp(1i * pi / 6)];
% alpha = [1; 0.9*exp(1i * pi / 4); 0.7*exp(1i * pi / 6)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars, alpha);
A = DOA.Array_Manifold(sensor_locations_oca, doa);
SNR_dB = -10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

%% Start The Algorithm 

% Ry = (1 / L) * (y * y');
Rs = alpha * alpha';
SNR = 10^(SNR_dB / 10);
Ry = A * Rs * A'% + (1 / SNR) * eye(M);

r = Ry(:);
z = Rearrange_According_to_Sensor_Locations(DOA, r, sensor_locations_oca);

R = zeros(N);
for i = 1:N
    start_idx = N + 1 - i;
    end_idx = start_idx + N - 1;
    R(:, i) = z(start_idx:end_idx);
end

holes = [14 16 18 19];
W = ones(N);
for i = 1:N
    for j = 1:N
        diff = abs(i - j);
        if ismember(diff, holes)
            W(i, j) = 0;
        end
    end
end

% Steepest Descent (does not work)

d = K; % # of columns in A and # of rows in B
A0 = 0.5 * (1 + 1i) * ones(N, d);
B0 = 0.5 * (1 + 1i) * ones(d, N);

LRMD = LRMD_Object(A0, B0, W, R, 0.5, 1);
[A_k1, B_k1] = LRMD.Steepest_Descent(100);

X = A_k1 * B_k1;

%%

% The optimization process starts here

beta = 1;
gamma = 0.5;
d = K; % # of columns in A and # of rows in B
N_max = 100;

A_k1 = 5 + 1i * 5 + 1 * (rand(N, d) + 1i * rand(N, d));
B_k1 = 5 + 1i * 5 + 1 * (rand(d, N) + 1i * rand(d, N));

lambda_A = 0.0001;
lambda_B = 0.0001;

for idx = 1:N_max
    A_k = A_k1;
    B_k = B_k1;
    % Update A
    for i = 1:N
        for j = 1:d
            grad = -0.5 * B_k(j, :) * (W(i, :) .* (R(i, :) - A_k(i, :) * B_k))';
            c = (A_k(i, j) - lambda_A * grad);
            A_k1(i, j) = (c / abs(c)) * max(abs(c) - gamma * lambda_A, 0);
        end
    end
    % Update B
    for i = 1:d
        for j = 1:N
            grad = -0.5 * (W(:, j) .* (R(:, j) - A_k1 * B_k(:, j)))' * A_k1(:, i);
            B_k1(i, j) = (1 / (1 + gamma * beta * lambda_B)) * (B_k(i, j) - lambda_B * grad);
        end
    end
end

X = A_k1 * B_k1; 

%% Plot The Spatial Spectrum

figure; 
ss = DOA.MUSIC(K, X, 0:N-1, angle_spec);
plot(angle_spec, 10*log10(ss));
xlabel('angle range (deg)')
ylabel('amplitude')
title('Spatial Spectrum')

% Elementary sonrasında Matris olarak işleme al

%% Plot The Spatial Spectrum

figure; 
ss = DOA.MUSIC(K, R(1:14, 1:14), 0:13, angle_spec);
plot(angle_spec, 10*log10(ss));
xlabel('angle range (deg)')
ylabel('amplitude')
title('Spatial Spectrum')

%% Succeeded

N = 5; d = 3;
A = (1 / sqrt(2)) * (rand(N, d) + 1i * rand(N, d));
D = diag(rand(d, 1)); 
B = D * A';

X = A * B;
rank(X)

W = ones(N);
W(4, 1) = 0;
W(5, 2) = 0;
W(1, 4) = 0;
W(2, 5) = 0;

R = W .* X;

beta = 1;
gamma = 0.5;
N_max = 100;

A_k1 = A + ((0.1 * rand(N, d) - 0.05) + 1i * (0.1 * rand(N, d) - 0.05));
B_k1 = B + ((0.1 * rand(d, N) - 0.05) + 1i * (0.1 * rand(d, N) - 0.05));

lambda_A = 0.0001;
lambda_B = 0.0001;

for idx = 1:N_max
    A_k = A_k1;
    B_k = B_k1;
    % Update A
    for i = 1:N
        for j = 1:d
            grad = -0.5 * B_k(j, :) * (W(i, :) .* (R(i, :) - A_k(i, :) * B_k))';
            c = (A_k(i, j) - lambda_A * grad);
            A_k1(i, j) = (c / abs(c)) * max(abs(c) - gamma * lambda_A, 0);
        end
    end
    % Update B
    for i = 1:d
        for j = 1:N
            grad = -0.5 * (W(:, j) .* (R(:, j) - A_k1 * B_k(:, j)))' * A_k1(:, i);
            B_k1(i, j) = (1 / (1 + gamma * beta * lambda_B)) * (B_k(i, j) - lambda_B * grad);
        end
    end
end

X = A_k1 * B_k1
R

%% Incoherent Case

N = 22; d = 2;
A = DOA.Array_Manifold(0:N-1, doa);
Rs = diag([1; 0.8]);
B = Rs * A';

X = A * B;
rank(X)
X_og = X;

holes = [14 16 18 19];
W = ones(N);
for i = 1:N
    for j = 1:N
        diff = abs(i - j);
        if ismember(diff, holes)
            W(i, j) = 0;
        end
    end
end

R = W .* X;

beta = 1;
gamma = 0.5;
N_max = 100;

A_k1 = A + ((0.1 * rand(N, d) - 0.05) + 1i * (0.1 * rand(N, d) - 0.05));
B_k1 = B + ((0.1 * rand(d, N) - 0.05) + 1i * (0.1 * rand(d, N) - 0.05));

A_0 = A_k1;
B_0 = B_k1;

% A_k1 = (2 * rand(N, d) - 1) + 1i * (2 * rand(N, d) - 1);
% B_k1 = (2 * rand(d, N) - 1) + 1i * (2 * rand(d, N) - 1);

lambda_A = 0.0001;
lambda_B = 0.0001;

for idx = 1:N_max
    A_k = A_k1;
    B_k = B_k1;
    % Update A
    for i = 1:N
        for j = 1:d
            grad = -0.5 * B_k(j, :) * (W(i, :) .* (R(i, :) - A_k(i, :) * B_k))';
            c = (A_k(i, j) - lambda_A * grad);
            A_k1(i, j) = (c / abs(c)) * max(abs(c) - gamma * lambda_A, 0);
        end
    end
    % Update B
    for i = 1:d
        for j = 1:N
            grad = -0.5 * (W(:, j) .* (R(:, j) - A_k1 * B_k(:, j)))' * A_k1(:, i);
            B_k1(i, j) = (1 / (1 + gamma * beta * lambda_B)) * (B_k(i, j) - lambda_B * grad);
        end
    end
end

X = A_k1 * B_k1
R

%%

% alpha = [1; (0.75+0.25*rand(K_coherent-1, 1)) .* exp(1i * pi * rand(K_coherent-1, 1))];
% % alpha = eye(K_coherent);
% B = DOA.kronecker(conj(A), A);
% P = alpha * alpha';
% p = P(:);
% 
% Bp = B * p;
% PHI = Rearrange_According_to_Sensor_Locations(DOA, Bp, sensor_locations_oca);
% 
% R = zeros(N);
% for i = 1:N
%     start_idx = N + 1 - i;
%     end_idx = start_idx + N - 1;
%     R(:, i) = PHI(start_idx:end_idx);
% end
% 
% %%
% 
% B22 = Rearrange_According_to_Sensor_Locations(DOA, B, sensor_locations_oca);
% 
% %%
% 
% R_temp = reshape(Bp, M, M);
% R2 = toeplitz(Virtual_Covariance_Column(DOA, R_temp, sensor_locations_oca)');
% 
% A2 = DOA.Array_Manifold(0:N-1, doa);
% % A2([15 17 19 20], :) = 0;
% % R_out = R_Toeplitz(A2 * P * A2', "full");
% % R_out = A2 * P * A2';
% R_out = A2 * diag(diag(P)) * A2';
% 
% %%
% 
% R_out_2 = A2 * P * A2';
% R_out_2 = R_Toeplitz(R_out_2, "full");
% % R_out_2 = R_out_2 / R_out_2(1, 1);
% % R = R / R(1,1);
% % 
% % norm(R-R_out_2, "fro")
% 
% %%
% 
% uDOF = DOA.Uniform_Degrees_Of_Freedom(sensor_locations_oca);
% N2 = (uDOF+1)*0.5;
% A3 = DOA.Array_Manifold(0:N2-1, doa);
% 
% Rs_3 = ((A3' * A3) \ A3') * R(1:N2, 1:N2) * (A3 / (A3' * A3));
% R_out_3 = A2 * diag(diag(Rs_3)) * A2';
% 
% %%
% 
% v3 = R_out_3(:, 1);
% v3([15 17 19 20]) = 0;
% R_out_3 = toeplitz(v3');
% 
% %%
% % R_out_3 = R_out_3 / R_out_3(1, 1);
% % R = R / R(1,1);
% 
% norm(R / R(1,1) - R_out_3 / R_out_3(1, 1), "fro")
% 
% %%
% 
% virtual_sensor_locations = 0:N-1;
% virtual_sensor_locations([15 17 19 20]) = 0;
% virtual_sensor_locations = sort(virtual_sensor_locations, "ascend");
% virtual_sensor_locations = virtual_sensor_locations(5:end);
% 
% A4 = DOA.Array_Manifold(virtual_sensor_locations, doa);
% v4 = R(:, 1);
% v4 = v4(virtual_sensor_locations + 1);
% R_in = toeplitz(v4');
% 
% 
% Rs_4 = ((A4' * A4) \ A4') * R_in * (A4 / (A4' * A4));
% R_out_4 = A2 * diag(diag(Rs_4)) * A2';
% 
% %%
% 
% figure; hold on
% ss = DOA.MUSIC(K, R(1:N2, 1:N2), 0:N2-1, angle_spec);
% plot(angle_spec, 10*log10(ss))
% ss = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
% plot(angle_spec, 10*log10(ss))
% ss = DOA.MUSIC(K, R_out_2, 0:N-1, angle_spec);
% plot(angle_spec, 10*log10(ss))
% ss = DOA.MUSIC(K, R_out_3, 0:N-1, angle_spec);
% plot(angle_spec, 10*log10(ss))
% ss = DOA.MUSIC(K, R_out_4, 0:N-1, angle_spec);
% plot(angle_spec, 10*log10(ss))
% legend('R', 'R_out', 'R_out_2', 'R_out_3', 'R_out_4')
% % legend('R_out', 'R_out_3')
% 
% %%
% 
% PHI3 = Rearrange_According_to_Sensor_Locations2(DOA, Bp, sensor_locations_oca);
% 
% R3 = zeros(N);
% for i = 1:N
%     start_idx = N + 1 - i;
%     end_idx = start_idx + N - 1;
%     R3(:, i) = PHI3(start_idx:end_idx);
% end

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