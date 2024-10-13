%% Notes for Future

% For normalization, you can use mean and variance of the vectorized form
% of the covariance matrix

%% Initialization

clear; clc; close all;

load Dataset.mat
pre_features = features;
clear features

%% Sensor Properties

sensor_locations = [0 1 4 7 9]; % SLA with 5 sensors

%% Dataset Preparation CRN2 (2 Source)

M = length(sensor_locations);
N = sensor_locations(M) + 1;

numOfData = size(labels, 2);

features = zeros(N, N, 3, numOfData);

for idx = 1:numOfData
    feature = pre_features(:, idx);
    feature = feature(1:N) + [0; feature(N+1:end)];
    R = toeplitz([feature(1) feature(2:end)']);
    features(:, :, 1, idx) = real(R);
    features(:, :, 2, idx) = imag(R);
    features(:, :, 3, idx) = angle(R) / pi;
end