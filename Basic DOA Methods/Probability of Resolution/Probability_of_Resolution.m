%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
sensor_locations = [0 1 4 10 16 18 21 23]; % MRA with 8 sensors
M = length(sensor_locations);
K = 2;      % # of sources
L = 100;    % # of snapshots
SNR_dB = 10;

phi_min = 30;
phi_max = 150;
delta_phi = 1;
doa = zeros(1, K);

noOfMethods = 4;
EPOCHS = 1000;
angle_seperation = 0.1:0.1:5;
prob_of_res = zeros(noOfMethods, length(angle_seperation));

angle_spec = phi_min:0.05:phi_max;
RMSE_Higher_Limit = 0.5; % const
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

%%

abs(A(:, 2)' * [A DOA.Array_Manifold(0.5, sensor_locations, 80)])

%%

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

K = 2;
L = 1000;
s = zeros(K, L);
freq_comp = 1;
for i = 1:freq_comp
    s = s + sqrt(1 / freq_comp) * DOA.Source_Generate(K, L);
end
(1 / L) * (s * s')

%% Plot The Graph

figure; hold on; grid on;
plot(angle_seperation, prob_of_res(1, :), 'b--o');
plot(angle_seperation, prob_of_res(2, :), 'r*');
plot(angle_seperation, prob_of_res(3, :));
plot(angle_seperation, prob_of_res(4, :));
xlabel("Angle Seperation"); ylabel("Probability of Resolution");
legend('CBF', 'Capon', 'MUSIC', 'DML');
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