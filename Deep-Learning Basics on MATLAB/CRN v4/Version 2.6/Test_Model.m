%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2_6_K3.mat

%% Test Case I

Test(net, DOA)

%%

function Test(net, DOA)
sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);

feature = zeros(M, M, 3);

delta_phi = 1;
phi_min = 30;
phi_max = 150;
K = 3;
doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

A_ohm = DOA.Array_Manifold(0.5, sensor_locations, doa);
L = 70; % # of snapshots
s = DOA.Source_Generate(K, L);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;
R_ohm = (1 / L) * (y * y');
normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

feature(:, :, 1) = real(normalized_R_ohm);
feature(:, :, 2) = imag(normalized_R_ohm);
feature(:, :, 3) = angle(R_ohm) / pi;

spatial_spectrum_1 = predict(net, feature).';
spatial_spectrum_1 = spatial_spectrum_1 / max(spatial_spectrum_1);
angles = phi_min:delta_phi:phi_max;
spatial_spectrum_2 = DOA.MUSIC(K, 0.5, R_ohm, sensor_locations, angles);

figure; hold on;
plot(angles, spatial_spectrum_1);
plot(angles, spatial_spectrum_2);
legend('Net', 'MUSIC');
title_text = "Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end