clear; clc; close all;

DOA = FunctionsOfDOA();
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
doa = 60;     % source angles
n = length(doa);    % # of sources
N = 100;            % # of snapshots
vars = 1;       % variance values of the source signals
SNR_dB_vals = -20:1:20;
EPOCHS = 1000;

noOfMethods = 4;
RMSE = zeros(noOfMethods, length(SNR_dB_vals));

for epoch = 1:EPOCHS
    s = DOA.Source_Generate(n, N, vars);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    C = DOA.Mutual_Coupling(100, 0.1, M, sensor_locations);
    for idx = length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        v = DOA.Noise_Generate(SNR_dB, M, N);
        y = C * A * s + v;

        % CBF
        doa_est = CBF(doa-5:0.01:doa+5, sensor_locations, y);
        RMSE(1, idx) = RMSE(1, idx) + (1 / EPOCHS) * (doa_est - doa)^2;

        Ry = (1 / N) * (y * y');

        % Capon
        doa_est = Capon(doa-5:0.01:doa+5, sensor_locations, y, Ry);
        RMSE(2, idx) = RMSE(2, idx) + (1 / EPOCHS) * (doa_est - doa)^2;

        % MUSIC
        doa_est = MUSIC(doa-5:0.01:doa+5, sensor_locations, Ry, M, n);
        RMSE(3, idx) = RMSE(3, idx) + (1 / EPOCHS) * (doa_est - doa)^2;

        % DML
        doa_est = DML(doa-5:0.01:doa+5, sensor_locations, Ry, M);
        RMSE(4, idx) = RMSE(4, idx) + (1 / EPOCHS) * (doa_est - doa)^2;
    end
end
RMSE = sqrt(RMSE);

figure; hold on;
plot(SNR_dB_vals, RMSE(1, :), 'b--o');
plot(SNR_dB_vals, RMSE(2, :), 'r*');
plot(SNR_dB_vals, RMSE(3, :));
plot(SNR_dB_vals, RMSE(4, :));
xlabel("SNR (dB)"); ylabel("RMSE");
legend('CBF', 'Capon', 'MUSIC', 'DML');

% CBF
function doa_est = CBF(theta, sensor_locations, y)
spec = zeros(length(theta), 1);
for i = 1:length(theta)
    h = exp(1i * pi * sensor_locations.' * cosd(theta(i)));
    spec(i) = abs(h' * (y * y') * h);
end
[~, max_idx] = max(spec);
doa_est = theta(max_idx);
end

% Capon
function doa_est = Capon(theta, sensor_locations, y, Ry)
spec = zeros(length(theta), 1);
for i = 1:length(theta)
    a = exp(1i * pi * sensor_locations.' * cosd(theta(i)));
    h = (Ry \ a) / (a' / Ry * a);
    spec(i) = abs(h' * (y * y') * h);
end
[~, max_idx] = max(spec);
doa_est = theta(max_idx);
end

% MUSIC
function doa_est = MUSIC(theta, sensor_locations, Ry, M, n)
[eig_vecs, ~] = eig(Ry);
U_N = eig_vecs(:, 1:M-n);
spec = zeros(length(theta), 1);
for i = 1:length(theta)
    a = exp(1i * pi * sensor_locations.' * cosd(theta(i)));
    spec(i) = 1 / abs(a' * (U_N * U_N') * a);
end
[~, max_idx] = max(spec);
doa_est = theta(max_idx);
end

% DML
function doa_est = DML(theta, sensor_locations, Ry, M)
spec = zeros(length(theta), 1);
for i = 1:length(theta)
    a = exp(1i * pi * sensor_locations.' * cosd(theta(i)));
    PI_a = a * (1 / (a' * a) * a');
    PI_a_ort = eye(M) - PI_a;

    spec(i) = trace(PI_a_ort * Ry);
end
[~, min_idx] = min(spec);
doa_est = theta(min_idx);
end