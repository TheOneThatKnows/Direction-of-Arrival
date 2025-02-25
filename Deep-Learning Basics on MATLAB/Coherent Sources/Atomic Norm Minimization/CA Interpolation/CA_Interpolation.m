%% Article Name: 
% DOA Estimation of Coherent Signals via Coprime Array Interpolation 
% https://ieeexplore.ieee.org/abstract/document/9044735

%% Initialize

clear; clc; close all;
addpath('E:\Direction-of-Arrival');
DOA = FunctionsOfDOA();

%% Define The Sensor Array

M = 3; N = 5;
sensor_locations_1 = 0:M:(N-1)*M;
sensor_locations_2 = 0:N:(M-1)*N;

sensor_locations = sort([sensor_locations_1 sensor_locations_2(2:end)]);

M = length(sensor_locations);
N = sensor_locations(M) + 1;

%% Simulate The Environment

K = 3;          % # of sources
K_coherent = 3;
L = 500;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

% doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
% doa = [50 120];
doa = [90 111.5 140];

vars = ones(K - K_coherent + 1, 1);
% vars(1:K_coherent) = [1; 0.75+0.25*rand(K_coherent-1, 1)];
% alpha = [1; exp(1i * pi / 6)];
alpha = [1; exp(1i * pi / 6); exp(1i * pi / 10)];
s = DOA.Source_Generate_Final(K, K_coherent, L, vars, alpha);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

y2 = zeros(N, L);
y2(sensor_locations+1, :) = y;

Ryu = (1 / L) * (y2 * y2');
m = 1;
N2 = 0.5 * (N + 1);
F_hat = zeros(N2);
C = zeros(N2);
for i = 1:N2
    for j = 1:N2
        diff = i - j;
        ind = N2 - diff;
        F_hat(i, j) = Ryu(m, ind);

        if ismember(ind-1, sensor_locations)
            C(i, j) = 1;
        end
    end
end

% does not work
F_tilda = solve_optimization(C, F_hat);

%% Plot The Spatial Spectrum

figure; 
ss = DOA.MUSIC(K, F_tilda, 0:N2-1, angle_spec);
plot(angle_spec, 10*log10(ss));
xlabel('angle range (deg)')
ylabel('amplitude')
title('Spatial Spectrum')