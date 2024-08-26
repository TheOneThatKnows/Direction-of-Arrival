%% Initialization

clear; clc; close all;

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

load CRN_Network_v2_3_v2.mat
crn_net_v2_3 = net;

clear net

%% Test Case I

Test(crn_net_v2_3, DOA)

%%

function Test(net, DOA)
sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(M) + 1;

feature = zeros(M, M, 2);

delta_phi = 3;
phi_min = 30;
phi_max = 150;
K = 2;
while true
    doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
    if (doa(2) - doa(1) > 2 * delta_phi)
        break
    end
end

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

angles_2 = phi_min:delta_phi:phi_max;
spec_2 = predict(net, feature(:, :, 1), feature(:, :, 2)).';

angles_1 = 30:1:150;
spec_1 = DOA.MUSIC(K, 0.5, R_ohm, sensor_locations, angles_1);

figure; hold on;
% plot(angles_1, spec_1);
% stem(angles_2, spec_2 / max(spec_2));
stem(angles_2, spec_2);
% legend('MUSIC', 'CRN v2.3');
title_text = "MUSIC Spectrum Comparison, DOAs: ";
for i = 1:K
    title_text = title_text + " " + doa(i);
end
title(title_text);
xlabel('angles (deg)');
ylabel('Amplitude');

DOA_Estimator(spec_2, angles_2)
end

function doa_est = DOA_Estimator(spec, angle_spec)
spec = spec(:); spec = [0; spec; 0];
angle_spec = angle_spec(:); angle_spec = [0; angle_spec; 0];
doa_est = zeros(1, 2);
[mags, inds] = findpeaks(spec);
[~, ind] = max(mags);
temp_idx1 = inds(ind);

mags = [mags(1:ind-1); mags(ind+1:end)];
inds = [inds(1:ind-1); inds(ind+1:end)];
[~, ind] = max(mags);
temp_idx2 = inds(ind);

idx1 = min(temp_idx1, temp_idx2);
idx2 = max(temp_idx1, temp_idx2);

if idx2 - idx1 == 2
    sum_first_half = sum(spec(idx1-1:idx1));
    sum_second_half = sum(spec(idx2:idx2+1));

    ratio1 = (100 - sum_first_half) / (200 - sum_first_half - sum_second_half);
    ratio2 = 1 - ratio1;

    doa_est(1) = [spec(idx1-1:idx1).' (spec(idx1+1) * ratio1)] * angle_spec(idx1-1:idx1+1) / sum([spec(idx1-1:idx1); (spec(idx1+1) * ratio1)]);
    doa_est(2) = [spec(idx2:idx2+1).' (spec(idx2-1) * ratio2)] * angle_spec(idx2-1:idx2+1) / sum([spec(idx2:idx2+1); (spec(idx2-1) * ratio2)]);

    return
end

doa_est(1) = spec(idx1-1:idx1+1).' * angle_spec(idx1-1:idx1+1) / sum(spec(idx1-1:idx1+1));
doa_est(2) = spec(idx2-1:idx2+1).' * angle_spec(idx2-1:idx2+1) / sum(spec(idx2-1:idx2+1));

end

function doa_est = DOA_Estimator2(spec, angle_spec)
spec = spec(:);
angle_spec = angle_spec(:);
doa_est = zeros(1, 2);
[mags, inds] = findpeaks(spec);
[~, ind] = max(mags);
idx = inds(ind);

if length(mags) == 1
    idx_init = max(1, idx - 2);
    idx_end = min(length(spec), idx + 2);
    sum_first_half = sum(spec(idx_init:idx-1));
    sum_second_half = sum(spec(idx+1:idx_end));
    ratio1 = (100 - sum_first_half) / (200 - sum_first_half - sum_second_half);
    ratio2 = 1 - ratio1;
    doa_est(1) = (spec(idx_init:idx-1).' * angle_spec(idx_init:idx-1) + mags * ratio1 * angle_spec(idx)) / (sum([spec(idx_init:idx-1); mags * ratio1]));
    doa_est(2) = (spec(idx+1:idx_end).' * angle_spec(idx+1:idx_end) + mags * ratio2 * angle_spec(idx)) / (sum([spec(idx+1:idx_end); mags * ratio2]));

    return
end

idx1 = idx;
mags = [mags(1:ind-1); mags(ind+1:end)];
inds = [inds(1:ind-1); inds(ind+1:end)];
[~, ind] = max(mags);
idx2 = inds(ind);

temp1 = idx1;
temp2 = idx2;
idx1 = min(temp1, temp2);
idx2 = max(temp1, temp2);

if idx2 - idx1 == 1


    return
end

if idx2 - idx1 == 2

    return
end
end