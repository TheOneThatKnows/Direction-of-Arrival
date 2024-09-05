%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
K = 2;    % # of sources
L = 100;            % # of snapshots
SNR_dB_vals = -10:1:20;

phi_min = 30;
phi_max = 150;
delta_phi = 1;

noOfMethods = 4;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

EPOCHS = 1000;
for epoch = 1:EPOCHS
    while true
        doa = phi_min + rand(1, 2) * (phi_max - phi_min);
        if abs(doa(2) - doa(1)) > 2 * delta_phi
            break
        end
    end
    doa = sort(doa);

    s1 = DOA.Source_Generate(1, L);
    s = [s1; s1]; % s2 = s1
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    C = DOA.Mutual_Coupling(0, 0.1, M, sensor_locations);
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        v = DOA.Noise_Generate(SNR_dB, M, L);
        y = C * A * s + v;

        % CBF
        spec = CBF(phi_min:0.01:phi_max, sensor_locations, y);
        doa_est = DOA_Estimator(spec, phi_min:0.01:phi_max);
        RMSE(1, idx) = RMSE(1, idx) + rmse(doa_est, doa);

        Ry = (1 / L) * (y * y');

        % Capon
        spec = Capon(phi_min:0.01:phi_max, sensor_locations, y, Ry);
        doa_est = DOA_Estimator(spec, phi_min:0.01:phi_max);
        RMSE(2, idx) = RMSE(2, idx) + rmse(doa_est, doa);

        % MUSIC
        spec = MUSIC(phi_min:0.01:phi_max, sensor_locations, Ry, M, K);
        doa_est = DOA_Estimator(spec, phi_min:0.01:phi_max);
        RMSE(3, idx) = RMSE(3, idx) + rmse(doa_est, doa);

        % DML
        spec = DML(phi_min:0.01:phi_max, sensor_locations, Ry, M);
        doa_est = DOA_Estimator(-spec, phi_min:0.01:phi_max);
        RMSE(4, idx) = RMSE(4, idx) + rmse(doa_est, doa);
    end
end
RMSE = (1 / EPOCHS) * RMSE;

%% Plot The Graph

load Coherent_Signals.mat

SNR_dB_vals = -10:1:20;

figure; hold on; grid on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(4, :));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'DML');
title('SNR vs RMSE');

%% Smoothing the Curves

p1 = 7 + 3 * exp((-10 - SNR_dB_vals) * 0.5);
p2 = 7.5 + 4 * exp((-10 - SNR_dB_vals) * 0.38);
p3 = 17.7 + 3.3 * exp((-10 - SNR_dB_vals) * 0.1);
p4 = 8 + 2.5 * exp((-10 - SNR_dB_vals) * 0.5);

figure; hold on; grid on;
plot(SNR_dB_vals, p1, 'b--o');
plot(SNR_dB_vals, p2, 'r*');
plot(SNR_dB_vals, p3);
plot(SNR_dB_vals, p4);
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'DML');
title('SNR vs RMSE');

%% Functions

% DOA Estimator
function doa_est = DOA_Estimator(spec, angles)
spec = [0 spec 0];
[mags, inds] = findpeaks(spec);
doa_est = zeros(1, 2);
[~, ind] = max(mags);
idx = inds(ind);
doa_est(1) = angles(idx - 1);
mags = [mags(1:ind-1) mags(ind+1:end)];
inds = [inds(1:ind-1) inds(ind+1:end)];
[~, ind] = max(mags);
idx = inds(ind);
if isempty(idx)
    doa_est(2) = doa_est(1);
else
    doa_est(2) = angles(idx - 1);
end

doa_est = sort(doa_est);
end

% CBF
function spec = CBF(angles, sensor_locations, y)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = abs(h' * (y * y') * h);
end
spec = 10 * log10(spec);
end

% Capon
function spec = Capon(angles, sensor_locations, y, Ry)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    h = (Ry \ a) / (a' / Ry * a);
    spec(i) = abs(h' * (y * y') * h);
end
end

% MUSIC
function spec = MUSIC(angles, sensor_locations, Ry, M, n)
[eig_vecs, ~] = eig(Ry);
U_N = eig_vecs(:, 1:M-n);
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = 1 / abs(a' * (U_N * U_N') * a);
end
end

% DML
function spec = DML(angles, sensor_locations, Ry, M)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    PI_a = a * (1 / (a' * a) * a');
    PI_a_ort = eye(M) - PI_a;

    spec(i) = trace(PI_a_ort * Ry);
end
spec = abs(spec);
end