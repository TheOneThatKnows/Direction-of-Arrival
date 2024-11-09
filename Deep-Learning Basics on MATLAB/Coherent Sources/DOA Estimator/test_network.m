%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2_6_K6.mat

%% Test Case I

[ss_crn, ss_music, angles] = Test(net, DOA);

%%

function [spatial_spectrum_1, spatial_spectrum_2, angles] = Test(net, DOA)
sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(end) + 1;

feature = zeros(N, N, 2);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
K = 5;
K_coherent = 2;
doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A_ohm = DOA.Array_Manifold(sensor_locations, doa);
L = 70; % # of snapshots
s = [DOA.Coherent_Source_Generate(K_coherent, L);
    DOA.Source_Generate(K - K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;
R_ohm = (1 / L) * (y * y');

z = R_ohm(:);
z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
R_out = zeros(N);
for i = 1:N
    z1_i = z1(i:i + N - 1);
    R_out = R_out + (1 / N) * (z1_i * z1_i');
end

R_norm = (R_out - mean(R_out(:))) / std(R_out(:));

feature(:, :, 1) = real(R_norm);
feature(:, :, 2) = imag(R_norm);

spatial_spectrum_1 = predict(net, feature).';
spatial_spectrum_1 = spatial_spectrum_1 / max(spatial_spectrum_1);

angles = phi_min:delta_phi:phi_max;
spatial_spectrum_2 = DOA.MUSIC(K, R_out, 0:N-1, angles);

figure; hold on;
plot(angles, spatial_spectrum_1);
plot(angles, spatial_spectrum_2);
legend('Net', 'SS-MUSIC');
title_text = "Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end