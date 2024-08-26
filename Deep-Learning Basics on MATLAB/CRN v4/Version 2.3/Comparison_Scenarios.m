%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2_2.mat

%% 
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

SNR_dB_vals = -10:1:10;
EPOCHS = 5000;

noOfMethods = 4;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angles = phi_min:delta_phi:phi_max;

for epoch = 1:EPOCHS
    while true
        doa = phi_min + rand(1, 2) * (phi_max - phi_min);
        if abs(doa(2) - doa(1)) > delta_phi
            break
        end
    end
    doa = sort(doa);

    s = DOA.Source_Generate(K, L);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y = A * s + n;

        % CBF
        spec = CBF(angles, sensor_locations, y);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(1, idx) = RMSE(1, idx) + rmse(doa_est, doa);

        Ry = (1 / L) * (y * y');

        % Capon
        spec = Capon(angles, sensor_locations, y, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(2, idx) = RMSE(2, idx) + rmse(doa_est, doa);

        % MUSIC
        spec = MUSIC(angles, sensor_locations, Ry, M, K);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(3, idx) = RMSE(3, idx) + rmse(doa_est, doa);

        % CRN_2 Network
        spec = CRN2_Function(net, M, Ry);
        doa_est = DOA_Estimator(spec, angles);
        doa_est = sort(doa_est);
        RMSE(4, idx) = RMSE(4, idx) + rmse(doa_est, doa);
    end
    if rem(epoch, 10) == 0
        disp(epoch)
    end
end

% RMSE = (1 / EPOCHS) * RMSE;
temp_RMSE = RMSE;
RMSE = (1 / (epoch - 1)) * RMSE;

figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(4, :));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'CRN_2 Network + MUSIC');
title('SNR vs RMSE')

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

function doa_est = DOA_Estimator_old(spec, angles, doa)
[mags, inds] = findpeaks(spec);
doa_est = zeros(1, 2);
[~, ind] = max(mags);
idx = inds(ind);
doa_est(1) = angles(idx);
mags = [mags(1:ind-1) mags(ind+1:end)];
inds = [inds(1:ind-1) inds(ind+1:end)];
[~, ind] = max(mags);
idx = inds(ind);
if isempty(idx)
    doa_est(2) = doa_est(1);
else
    doa_est(2) = angles(idx);
end

[~, ind] = min(abs(doa - doa_est(1)));
if ind == 2
    doa = doa(2:-1:1);
end

if rmse(doa, [doa_est(1) doa_est(1)]) < rmse(doa, doa_est)
    doa_est = [doa_est(1) doa_est(1)];
end
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
end

% CRN2

function spec = CRN2_Function(net, M, R_ohm)
normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

feature = zeros(M, M, 2);
feature(:, :, 1) = real(normalized_R_ohm);
feature(:, :, 2) = imag(normalized_R_ohm);

spec = predict(net, feature(:, :, 1), feature(:, :, 2));
spec = spec / max(spec);
end