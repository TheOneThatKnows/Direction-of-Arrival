%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load 'Comparison Dataset.mat'

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
EPOCHS = size(dataset_y, 4);

noOfMethods = 8;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angle_spec = phi_min:delta_phi:phi_max;

uDOF = DOA.Uniform_Degrees_Of_Freedom(sensor_locations);
M_v = 0.5 * (uDOF + 1);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angle_spec);
A_d = DOA.khatri_rao(conj(A_sparse), A_sparse);
[A1, ~] = DOA.Rearrange_According_to_Sensor_Locations(A_d, sensor_locations);
A2 = [A1(1:M_v-1, :); A1(M_v+1:end, :)];
I = eye(M); I = I(:);
I2 = zeros(uDOF, 1); I2(M_v) = 1;

for epoch = 1:EPOCHS
    doa = dataset_doa(:, epoch).';

    for idx = 1:length(SNR_dB_vals)
        y = dataset_y(:, :, idx, epoch);
        method = 1;

        % CBF
        spec = CBF(angle_spec, sensor_locations, y);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        Ry = (1 / L) * (y * y');

        % Capon
        spec = Capon(angle_spec, sensor_locations, y, Ry);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % MUSIC
        spec = MUSIC(angle_spec, sensor_locations, Ry, M, K);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % SS-MUSIC
        z = Ry(:);
        z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
        R_z1 = zeros(N);
        for i = 1:N
            z1_i = z1(i:i + N - 1);
            R_z1 = R_z1 + (1 / N) * (z1_i * z1_i');
        end
        spec = MUSIC(angle_spec, 0:N-1, R_z1, N, K);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % CS1
        spec = CS1(y, A_sparse, L).';
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % CS_Framework_Utilizing_Difference_Coarray 1
        z = Ry(:);
        spec = CS_Framework_Utilizing_Difference_Coarray_1(z, A_d, I).';
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % CS_Framework_Utilizing_Difference_Coarray 2
        [z1, ~] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
        spec = CS_Framework_Utilizing_Difference_Coarray_2(z1, A1, I2).';
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % CS_Framework_Utilizing_Difference_Coarray 3
        z2 = [z1(1:M_v-1); z1(M_v+1:end)];
        spec = CS_Framework_Utilizing_Difference_Coarray_3(z2, A2).';
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;
    end
    if rem(epoch, 10) == 0
        disp(epoch)
    end
end

RMSE = (1 / EPOCHS) * RMSE;
%%
figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(4, :));
plot(SNR_dB_vals, RMSE(5, :));
plot(SNR_dB_vals, RMSE(6, :));
plot(SNR_dB_vals, RMSE(7, :));
plot(SNR_dB_vals, RMSE(8, :));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'SS-MUSIC', 'CS_1', 'CS_{Sparse}_1', 'CS_{Sparse}_2', 'CS_{Sparse}_3');
title_text = "SNR vs RMSE (K=" + K + ")";
title(title_text)

%%

figure; hold on;
plot(SNR_dB_vals, RMSE_CBF, 'b--o');
plot(SNR_dB_vals, RMSE_Capon, 'r*');
plot(SNR_dB_vals, RMSE_MUSIC);
plot(SNR_dB_vals, RMSE_SS_MUSIC);
plot(SNR_dB_vals, RMSE_CS1);
plot(SNR_dB_vals, RMSE_CS_Coarray_1);
plot(SNR_dB_vals, RMSE_CS_Coarray_2);
plot(SNR_dB_vals, RMSE_CS_Coarray_3);
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'SS-MUSIC', 'CS_1', 'CS_{Sparse}_1', 'CS_{Sparse}_2', 'CS_{Sparse}_3');
title_text = "SNR vs RMSE (K=" + K + ")";
title(title_text)

%% Functions

% DOA Estimator
function doa_est = DOA_Estimator(spec, angles, K)
spec = [0 spec 0];
[mags, inds] = findpeaks(spec);
doa_est = zeros(1, K);
[~, ind] = max(mags);
idx = inds(ind);
doa_est(1) = angles(idx - 1);

for i = 2:K
    mags = [mags(1:ind-1) mags(ind+1:end)];
    inds = [inds(1:ind-1) inds(ind+1:end)];
    [~, ind] = max(mags);
    idx = inds(ind);
    if isempty(idx)
        doa_est(i:K) = doa_est(i-1);
        break
    else
        doa_est(i) = angles(idx - 1);
    end
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

% CS_Framework

function spec = CS1(y, A_sparse, snapshots)
CS = CS_Framework(y, A_sparse);
CS = CS.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
x = var(CS.Sg.').';
x = abs(x);
x = x / max(x);
spec = 10*log10(x);
end

function spec = CS_Framework_Utilizing_Difference_Coarray_1(z, A_d, I)
CS = CS_Framework_Utilizing_Difference_Coarray(z, [A_d I], [1 1]);
CS = CS.Steepest_Descent_DOA();
x = abs(CS.sg(1:end-1));
x = x / max(x);
spec = 10 * log10(x);
end

function spec = CS_Framework_Utilizing_Difference_Coarray_2(z1, A1, I2)
CS = CS_Framework_Utilizing_Difference_Coarray(z1, [A1 I2], [1 1]);
CS = CS.Steepest_Descent_DOA();
x = abs(CS.sg(1:end-1));
x = x / max(x);
spec = 10 * log10(x);
end

function spec = CS_Framework_Utilizing_Difference_Coarray_3(z2, A2)
CS = CS_Framework_Utilizing_Difference_Coarray(z2, A2, [1 1]);
CS = CS.Steepest_Descent_DOA();
x = abs(CS.sg);
x = x / max(x);
spec = 10 * log10(x);
end