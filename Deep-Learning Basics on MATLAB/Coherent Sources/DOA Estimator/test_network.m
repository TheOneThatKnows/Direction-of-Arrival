%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Coherent Sources']);

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival' ...
    '\Deep-Learning Basics on MATLAB\Custom Layers']);

load Network_pt2.mat

%% Test Case I

[ss_crn, ss_music, angles] = Test(net, DOA);
%%
figure; hold on;
plot(angles, ss_crn);
plot(angles, ss_music);
%%

function [spatial_spectrum_1, spatial_spectrum_2, angles] = Test(net, DOA)
sensor_locations = 0:10;
M = length(sensor_locations);
N = sensor_locations(end) + 1;

feature = zeros(M, M, 3);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
K = 3;
K_coherent = 2;
doa = DOA.DOA_Generate(K, phi_min, phi_max, 2 * delta_phi);

A = DOA.Array_Manifold(sensor_locations, doa);
L = 70; % # of snapshots
s = [DOA.Coherent_Source_Generate(K_coherent, L);
    DOA.Source_Generate(K - K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;
Ry = (1 / L) * (y * y');

R_toeplitz = R_Toeplitz(Ry, "half");

R_norm = (R_toeplitz - mean(R_toeplitz(:))) / std(R_toeplitz(:));

feature(:, :, 1) = real(R_norm);
feature(:, :, 2) = imag(R_norm);
feature(:, :, 3) = angle(R_norm) / pi;

spatial_spectrum_1 = predict(net, feature).';
spatial_spectrum_1 = spatial_spectrum_1 / max(spatial_spectrum_1);

angles = phi_min:delta_phi:phi_max;
spatial_spectrum_2 = DOA.MUSIC(K, R_toeplitz, 0:N-1, angles);

figure; hold on;
plot(angles, 10*log10(spatial_spectrum_1));
plot(angles, 10*log10(spatial_spectrum_2));
legend('Net', 'MUSIC');
title_text = "Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end