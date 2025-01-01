%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\' ...
    'Deep-Learning Basics on MATLAB\Coherent Sources']);

load Network_pt2.mat

%% Test Case I

[ss_crn, ss_music, angles] = Test(net, DOA);
%%
figure; hold on;
plot(angles, ss_crn);
plot(angles, ss_music);
%%

function [spatial_spectrum_1, spatial_spectrum_2, angle_spec] = Test(net, DOA)
M = 4; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations_ca = sort([sensor_locations_1 sensor_locations_2(2:end)]);
sensor_locations_oca = [0 max(M,N)+sensor_locations_ca(2:end)];

M = length(sensor_locations_oca);
N = sensor_locations_oca(M) + 1;

% virtual sensors

virtual_sensor_locations = 0:N-1;
virtual_sensor_locations([15 17 19 20]) = 0;
virtual_sensor_locations = sort(virtual_sensor_locations, "ascend");
virtual_sensor_locations = virtual_sensor_locations(5:end);

feature = zeros(N, N, 2);

L = 70;        % # of snapshots
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution
Q = (phi_max - phi_min) / delta_phi + 1;

K = 5;
K_coherent = randi(K);
K_coherent = K_coherent * sign(K_coherent - 1)

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

% Array Manifold
A = DOA.Array_Manifold(sensor_locations_oca, doa);

% signal generate
vars = ones(K-(K_coherent-1)*sign(K_coherent), 1);
s = DOA.Source_Generate_Final(K, K_coherent, L, vars);

% noise generate
SNR_dB = -5 + sign(randi(3) - 2.5) * 5 + rand * 10;
n = DOA.Noise_Generate(SNR_dB, M, L);

% measurements
y = A * s + n;
R = (1 / L) * (y * y');

R_toeplitz = toeplitz(Virtual_Covariance_Column(DOA, R, sensor_locations_oca)');

R_norm = (R_toeplitz - mean(R_toeplitz(:))) / std(R_toeplitz(:));

feature(:, :, 1) = real(R_norm);
feature(:, :, 2) = imag(R_norm);

column1 = predict(net, feature).';
column1 = column1(:);
column1 = column1(1:N) + 1i * [0; column1(N+1:end)];
R_out = toeplitz(column1');

v = R_toeplitz(:, 1);
v = v(virtual_sensor_locations + 1);
R_in = toeplitz(v');

angle_spec = phi_min:delta_phi:phi_max;
spatial_spectrum_1 = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
spatial_spectrum_2 = DOA.MUSIC(K, R_in, virtual_sensor_locations, angle_spec);

figure; hold on;
plot(angle_spec, 10*log10(spatial_spectrum_1));
plot(angle_spec, 10*log10(spatial_spectrum_2));
legend('Net+MUSIC', 'MUSIC');
title_text = "Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude (dB)');
end