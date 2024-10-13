%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load 'Comparison Dataset.mat' dataset_doa

%% 
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
K = 2;          % # of sources
L = 70;        % # of snapshots

SNR_dB_vals = -10:1:10;
EPOCHS = size(dataset_doa, 2);

CRLB = zeros(1, length(SNR_dB_vals));

for epoch = 1:EPOCHS
    doa = dataset_doa(:, epoch).';

    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        SNR = 10^(SNR_dB / 10);
        CRLB(idx) = CRLB(idx) + M / (2 * L * SNR * sind(doa(1)^2 * pi^2 * ((sensor_locations) * (sensor_locations).')));
    end
end
CRLB = CRLB / EPOCHS;

plot(SNR_dB_vals, CRLB);
%%
for epoch = 1:EPOCHS
    doa = dataset_doa(:, epoch).';

    s = DOA.Source_Generate(K, L);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    PI_A = A * (eye(K) / (A' * A)) * A';
    PI_A_Ort = eye(M) - PI_A;
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y = A * s + n;

        Rs = (1 / L) * (s * s');
        Ry = (1 / L) * (y * y');

        D = -1i * pi * (sensor_locations).' * sind(doa) .* exp(1i * pi * (sensor_locations).' * cosd(doa));
        crlb = var(n(1, :)) / (2 * L) * (eye(K) / real((D' * PI_A_Ort * D) .* (Rs * A' * (eye(M) / Ry) * A * Rs)'));
    end
    if rem(epoch, 10) == 0
        disp(epoch)
    end
    break
end

CRLB = (1 / EPOCHS) * CRLB;
%%
figure; hold on;
plot(SNR_dB_vals, CRLB);
xlabel("SNR (dB)"); ylabel("RMSE");
title("Cramer-Rao Lower Bound")

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