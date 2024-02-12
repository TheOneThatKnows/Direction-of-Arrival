%% Paper
% Name:
% Deep Neural Networks for Direction of Arrival Estimation of Multiple 
% Targets with Sparse Prior for Line-of-Sight Scenarios

% Link: https://ieeexplore.ieee.org/document/9963544


%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Start

% Sensor Array
sensor_locations = [0 1 2 3 7 11]; % ULA with M number of sensors
M = length(sensor_locations);

% Environment
doa = [86 97];
K = length(doa);
L = 500; % snapshots
SNR_dB = 10;
s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(coef, sensor_locations, doa);
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;
R = (1 / L) * (y * y');

possible_angles = linspace(80, 100, 210);
P = length(possible_angles);
A_star = zeros(M, P);
for p = 1:P
    A_star(:, p) = exp(1i * pi * (sensor_locations).' * cosd(possible_angles(p)));
end
A_tilda = DOA.khatri_rao(conj(A_star), A_star);

c = R(:);

ksa = A_tilda' * c;
figure;
plot(possible_angles, ksa);