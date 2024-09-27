%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival' ...
    '\Deep-Learning Basics on MATLAB\Custom Layers']);
addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival' ...
    '\Deep-Learning Basics on MATLAB\CRN v4\Version 2.6']);
load CRN_Network_v2_6.mat

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
% sensor_locations = [0 1 4 10 16 18 21 23]; % MRA with 8 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;      % # of sources
L = 70;    % # of snapshots
SNR_dB = 10;

phi_min = 30;
phi_max = 150;
delta_phi = 1;
doa = zeros(1, K);

noOfMethods = 6;
EPOCHS = 1000;
angle_seperation = 0.5:0.5:10;

%% Calculation

prob_of_res = zeros(noOfMethods, length(angle_seperation));

angle_spec = phi_min:1:phi_max;
RMSE_Higher_Limit = 2; % const
for i = 1:length(angle_seperation)
    for epoch = 1:EPOCHS
        % doa(1) = phi_min + rand * (phi_max - phi_min - angle_seperation(i));
        % doa(2) = doa(1) + angle_seperation(i);

        doa(1) = 90 - angle_seperation(i);
        doa(2) = 90 + angle_seperation(i);

        s = DOA.Source_Generate(K, L);
        A = DOA.Array_Manifold(0.5, sensor_locations, doa);
        C = DOA.Mutual_Coupling(0, 0.1, M, sensor_locations);
        v = DOA.Noise_Generate(SNR_dB, M, L);
        y = C * A * s + v;

        % CBF
        spec = CBF(angle_spec, sensor_locations, y);
        doa_est = DOA_Estimator(spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(1, i) = prob_of_res(1, i) + 1;
        end

        Ry = (1 / L) * (y * y');

        % Capon
        spec = Capon(angle_spec, sensor_locations, y, Ry);
        doa_est = DOA_Estimator(spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(2, i) = prob_of_res(2, i) + 1;
        end

        % MUSIC
        spec = MUSIC(angle_spec, sensor_locations, Ry, M, K);
        doa_est = DOA_Estimator(spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(3, i) = prob_of_res(3, i) + 1;
        end

        % DML
        spec = DML(angle_spec, sensor_locations, Ry, M);
        doa_est = DOA_Estimator(-spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(4, i) = prob_of_res(4, i) + 1;
        end

        % SS-MUSIC
        z = Ry(:);
        z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
        R_z1 = zeros(N);
        for j = 1:N
            z1_j = z1(j:j + N - 1);
            R_z1 = R_z1 + (1 / N) * (z1_j * z1_j');
        end
        spec = MUSIC(angle_spec, 0:N-1, R_z1, N, K);
        doa_est = DOA_Estimator(spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(5, i) = prob_of_res(5, i) + 1;
        end

        % CRN
        spec = CRN2_Function(net, M, Ry);
        doa_est = DOA_Estimator(spec, angle_spec);
        if rmse(doa_est, doa) < RMSE_Higher_Limit
            prob_of_res(6, i) = prob_of_res(6, i) + 1;
        end
    end
end

prob_of_res = (1 / EPOCHS) * prob_of_res;

%%

doa(1) = 90 - angle_seperation(end);
doa(2) = 90 + angle_seperation(end);

s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
C = DOA.Mutual_Coupling(0, 0.1, M, sensor_locations);
v = DOA.Noise_Generate(SNR_dB, M, L);
y = C * A * s + v;

% CBF
spec = CBF(angle_spec, sensor_locations, y);
doa_est = DOA_Estimator(spec, angle_spec);
plot(angle_spec, 10 * log10(spec/max(spec)))

%% Plot The Graph

figure; hold on; grid on;
plot(angle_seperation, prob_of_res(1, :), 'b--o');
plot(angle_seperation, prob_of_res(2, :), 'r--*');
plot(angle_seperation, prob_of_res(3, :), '--');
plot(angle_seperation, prob_of_res(4, :), '--');
plot(angle_seperation, prob_of_res(5, :), '--');
plot(angle_seperation, prob_of_res(6, :), '--');
xlabel("Angle Seperation"); ylabel("Probability of Resolution");
legend('CBF', 'Capon', 'MUSIC', 'DML', 'SS-MUSIC', 'DL Model');
title('Probability of Resolution');

%% Smoothing the Curves

spec = MUSIC(angle_spec, sensor_locations, Ry, M, K);
doa_est = DOA_Estimator(spec, angle_spec);

rmse(doa, doa_est)

figure;
plot(angle_spec, spec);

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

% CRN2

function spec = CRN2_Function(net, M, R_ohm)
normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

feature = zeros(M, M, 3);
feature(:, :, 1) = real(normalized_R_ohm);
feature(:, :, 2) = imag(normalized_R_ohm);
feature(:, :, 3) = angle(R_ohm) / pi;

spec = predict(net, feature);
spec = spec / max(spec);
end