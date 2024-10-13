%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load 'Comparison Dataset.mat'

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\CRN v4\Version 4.0']);
load CRN_Network_v4_0.mat
net_v4_0 = net;
clear net

load CRN_Network_v4_1.mat
net_v4_1 = net;
clear net

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
EPOCHS = size(dataset_doa, 2);

noOfMethods = 2;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

angle_spec = phi_min:delta_phi:phi_max;

A_sparse = DOA.Array_Manifold(sensor_locations, angle_spec);
A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);
feature_1 = zeros(M, M, 3);
for epoch = 1:EPOCHS
    doa = dataset_doa(:, epoch).';

    for idx = 1:length(SNR_dB_vals)
        y = dataset_y(:, :, idx, epoch);

        Ry = (1 / L) * (y * y');
        normalized_Ry = Ry / max(diag(abs(Ry)));
        feature_1(:, :, 1) = real(normalized_Ry);
        feature_1(:, :, 2) = imag(normalized_Ry);
        feature_1(:, :, 3) = angle(Ry) / pi;

        z = Ry(:);
        feature_2 = abs(A1' * z);

        method = 1;

        % v4.0
        spec = predict(net_v4_0, feature_1, feature_2.');
        spec = spec / max(spec);
        doa_est = DOA_Estimator(spec, angle_spec, K);
        RMSE(method, idx) = RMSE(method, idx) + rmse(doa_est, doa);
        method = method + 1;

        % v4.1
        spec = predict(net_v4_1, feature_1, feature_2.');
        spec = spec / max(spec);
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
for i = 1:noOfMethods
    plot(SNR_dB_vals, RMSE(i, :));
end
xlabel("SNR (dB)"); ylabel("RMSE");
legend('v4.0', 'v4.1');
title_text = "SNR vs RMSE (K=" + K + ")";
title(title_text)

%% Not used

figure; hold on;
plot(SNR_dB_vals, RMSE_CBF, 'b--o');
plot(SNR_dB_vals, RMSE_Capon, 'r*');
plot(SNR_dB_vals, RMSE_MUSIC);
plot(SNR_dB_vals, RMSE_SS_MUSIC);
plot(SNR_dB_vals, RMSE_DA_MUSIC, '*');
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'SS-MUSIC', 'DA-MUSIC');
title_text = "SNR vs RMSE (K=" + K + ")";
title(title_text)

%% Functions

% DOA Estimator
function doa_est = DOA_Estimator(spec, angles, K)
try
    spec = [-inf spec -inf];
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
catch
    doa_est = [90 90];
end
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
spec = spec / max(spec);
end

% Capon
function spec = Capon(angles, sensor_locations, y, Ry)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    h = (Ry \ a) / (a' / Ry * a);
    spec(i) = abs(h' * (y * y') * h);
end
spec = spec / max(spec);
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
spec = spec / max(spec);
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