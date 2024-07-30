%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2_2.mat
crn_net_v2_2 = net;

clear net

%% Test Case I

Test(crn_net_v2_2, DOA)

%%

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
L = 70; % # of snapshots
s = DOA.Source_Generate(K, L);
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A_ohm * s + n;
R_ohm = (1 / L) * (y * y');
normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

feature(:, :, 1) = real(normalized_R_ohm);
feature(:, :, 2) = imag(normalized_R_ohm);

spec_2 = predict(net, feature(:, :, 1), feature(:, :, 2)).';

angles = 30:1:150;
spec_1 = DOA.MUSIC(K, 0.5, R_ohm, sensor_locations, angles);

figure; hold on;
plot(angles, spec_1);
plot(angles, spec_2 / max(spec_2));
legend('MUSIC', 'CRN v2.2');
title_text = "MUSIC Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end