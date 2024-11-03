%% Article

% Title: 
% Direction-of-Arrival Estimation of Mixed Coherent and Uncorrelated Signals 

% Link:
% https://ieeexplore.ieee.org/abstract/document/10640208?casa_token=8V0l9DteRzYAAAAA:QzEjv5QNOgmTkZ3-I4KgSoOrpftMW9hfe6YO8UcnNKHG59qPupgeLd59_G4e9rG0Kj5Ig-fw6w

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 1)

s = [DOA.Source_Generate(K-K_coherent, L); 
    DOA.Coherent_Source_Generate(K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');

spec = DOA.MUSIC(K, R, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10 * log10(spec));

R_ = zeros(M, M, (M-1 - rem(M-1, 2))/2);
R_m = zeros(M);

C_ = zeros(M, M, (M-1 - rem(M-1, 2))/2);
C_m = zeros(M);

for m = 0:(M-1 - rem(M-1, 2))/2
    J = zeros(M - m);
    for j = 1:M-m
        J(j, M-m-j+1) = 1;
    end

    % R
    r_m = [zeros(1, m) (R(m+1:M, m+1).' * J) R(m+1, m+2:M) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        R_m(i+1, :) = r_m(start_idx:end_idx).';
    end

    R_(:, :, m+1) = R_m;

    % C
    c_m = [zeros(1, m) ones(1, 2*M-2*m-1) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        C_m(i+1, :) = c_m(start_idx:end_idx).';
    end

    C_(:, :, m+1) = C_m;
end

R_sum = zeros(M);
C_sum = zeros(M);
for m = 0:(M-1)/2
    R_sum = R_sum + R_(:, :, m+1);
    C_sum = C_sum + C_(:, :, m+1);
end

R_out = R_sum ./ C_sum;

spec = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
plot(angle_spec, 10 * log10(spec));
title_text = "DOAs: [";
for k = 1:K
    title_text = title_text + " " + doa(k);
end
title_text = title_text + " ]";
title(title_text)
legend('1', '2')

%% Initialization

clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

% sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
% sensor_locations = [1 2 3 4 8 12] - 1; % NA with 6 sensors
sensor_locations = [1 2 3 4 9 13] - 1; % NA v2 with 6 sensors
% sensor_locations = 0:10; % ULA with 11 sensors

M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 3;          % # of sources
K_coherent = 2;
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, 10)

s = [DOA.Source_Generate(K-K_coherent, L); 
    DOA.Coherent_Source_Generate(K_coherent, L)];
shuffledIndices = randperm(K);
s = s(shuffledIndices, :);
A = DOA.Array_Manifold(sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

R = (1 / L) * (y * y');

spec = DOA.MUSIC(K, R, sensor_locations, angle_spec);

figure; hold on
plot(angle_spec, 10 * log10(spec));

R = Spatial_Smoothing(DOA, sensor_locations, R);

M = size(R, 1);

R_ = zeros(M, M, (M-1 - rem(M-1, 2))/2);
R_m = zeros(M);

C_ = zeros(M, M, (M-1 - rem(M-1, 2))/2);
C_m = zeros(M);

for m = 0:(M-1 - rem(M-1, 2))/2
    J = zeros(M - m);
    for j = 1:M-m
        J(j, M-m-j+1) = 1;
    end

    % R
    r_m = [zeros(1, m) (R(m+1:M, m+1).' * J) R(m+1, m+2:M) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        R_m(i+1, :) = r_m(start_idx:end_idx).';
    end

    R_(:, :, m+1) = R_m;

    % C
    c_m = [zeros(1, m) ones(1, 2*M-2*m-1) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        C_m(i+1, :) = c_m(start_idx:end_idx).';
    end

    C_(:, :, m+1) = C_m;
end

R_sum = zeros(M);
C_sum = zeros(M);
for m = 0:(M-1)/2
    R_sum = R_sum + R_(:, :, m+1);
    C_sum = C_sum + C_(:, :, m+1);
end

R_out = R_sum ./ C_sum;

spec = DOA.MUSIC(K, R_out, 0:N-1, angle_spec);
plot(angle_spec, 10 * log10(spec));
title_text = "DOAs: [";
for k = 1:K
    title_text = title_text + " " + doa(k);
end
title_text = title_text + " ]";
title(title_text)
legend('1', '2')

%%

function R_z1 = Spatial_Smoothing(DOA, sensor_locations, R)
N = sensor_locations(end) - sensor_locations(1) + 1;
z = R(:);
z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
R_z1 = zeros(N);
for i = 1:N
    z1_i = z1(i:i + N - 1);
    R_z1 = R_z1 + (1 / N) * (z1_i * z1_i');
end
end