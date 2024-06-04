%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2.mat

%% Test Case I

Test(net, DOA)

%%

function R = r2R(r)
N = (length(r) + 1) * 0.5;
R_c1 = r(1:N) + 1i * [0; r(N+1:end)];
R = toeplitz(R_c1');
end

function Test(net, DOA)
sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(M) + 1;

feature = zeros(M, M, 2);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
K = 2;
doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A_ohm = DOA.Array_Manifold(0.5, sensor_locations, doa);
L = 100; % # of snapshots
s = DOA.Source_Generate(K, L);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;
R_ohm = (1 / L) * (y * y');
normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

feature(:, :, 1) = real(normalized_R_ohm);
feature(:, :, 2) = imag(normalized_R_ohm);

R_c1 = predict(net, feature(:, :, 1), feature(:, :, 2)).';
R = r2R(R_c1);

[spatial_spectrum_1, angles] = DOA.MUSIC(K, 0.5, R_ohm, sensor_locations);
[spatial_spectrum_2, ~] = DOA.MUSIC(K, 0.5, R, 0:N-1);

figure; hold on;
plot(angles, spatial_spectrum_1);
plot(angles, spatial_spectrum_2);
legend('R_ohm', 'R');
title_text = "MUSIC Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end