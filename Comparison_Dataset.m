%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

rng shuffle
%% 
sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

SNR_dB_vals = -10:1:10;
noOfTrials = 5000;

dataset_y = zeros(M, L, length(SNR_dB_vals), noOfTrials);
dataset_doa = zeros(K, noOfTrials);

for trial = 1:noOfTrials
    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

    s = DOA.Source_Generate(K, L);
    A = DOA.Array_Manifold(0.5, sensor_locations, doa);
    for idx = 1:length(SNR_dB_vals)
        SNR_dB = SNR_dB_vals(idx);
        n = DOA.Noise_Generate(SNR_dB, M, L);
        y = A * s + n;

        dataset_y(:, :, idx, trial) = y;
    end
    dataset_doa(:, trial) = doa.';
end