%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Encoder']);
load First_Level.mat
load Second_Level.mat

load Off_Grid_Network

%% Test Case I

[ss_crn, ss_music, angles] = Test(net, DOA);

%%

function [spatial_spectrum_1, spatial_spectrum_2, angles] = Test(net, DOA)
sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(end) + 1;

delta_phi = 1;
phi_min = 30;
phi_max = 150;
K = 2;
doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A_ohm = DOA.Array_Manifold(sensor_locations, doa);
L = 70; % # of snapshots
s = DOA.Source_Generate(K, L);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;
R_ohm = (1 / L) * (y * y');

[eig_vecs, ~] = eig(R_ohm);
noise_vecs = eig_vecs(:, 1:M-K);
G = noise_vecs * noise_vecs';

G_u = [real(G(:, 1)); imag(G(:, 1))];
mean_G = mean(G_u);

feature = (G_u - mean_G) / max(abs(G_u - mean_G));

spatial_spectrum_1 = predict(net, feature.').';
spatial_spectrum_1 = spatial_spectrum_1 / max(spatial_spectrum_1);

z = R_ohm(:);
z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
R_z1 = zeros(N);
for i = 1:N
    z1_i = z1(i:i + N - 1);
    R_z1 = R_z1 + (1 / N) * (z1_i * z1_i');
end
angles = phi_min:delta_phi:phi_max;
spatial_spectrum_2 = DOA.MUSIC(K, R_z1, 0:N-1, angles);

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