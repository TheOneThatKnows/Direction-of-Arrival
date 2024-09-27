%% CS_Framework

clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);

K = 2; L = 70; SNR_dB = 15;
phi_min = 30;
phi_max = 150;
delta_phi = 1;
angle_spec = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);
s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
C = DOA.Mutual_Coupling(0, 0.1, M, sensor_locations);
n = DOA.Noise_Generate(SNR_dB, M, L);
y = C * A * s + n;

Ry = (1 / L) * (y * y');
z = Ry(:);
[z1, ~] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);

A_sparse = DOA.Array_Manifold(0.5, sensor_locations, angle_spec);
A_d = DOA.khatri_rao(conj(A_sparse), A_sparse);
[A1, ~] = DOA.Rearrange_According_to_Sensor_Locations(A_d, sensor_locations);

I = eye(M); I = I(:);

uDOF = length(z1);
M_v = 0.5 * (uDOF + 1);
I2 = zeros(uDOF, 1); I2(M_v) = 1;

z2 = [z1(1:M_v-1); z1(M_v+1:end)];
A2 = [A1(1:M_v-1, :); A1(M_v+1:end, :)];

%% CS1

CS = CS_Framework(y, A_sparse);

CS = CS.Steepest_Descent_DOA([1 min(1/L, 0.1)]);

if L ~= 1
    x = var(CS.Sg.').';
    x = abs(x);
else
    x = abs(CS.Sg);
end
x = x / max(x);
figure; plot(angle_spec, 10*log10(x));
title_text = "DOA Angles: " + doa(1) + ", " + doa(2);
title(title_text)

%% CS_Framework_Utilizing_Difference_Coarray 1

CS = CS_Framework_Utilizing_Difference_Coarray(z, [A_d I], [1 1]);

CS = CS.Steepest_Descent_DOA();

x = abs(CS.sg(1:end-1));
x = x / max(x);
figure; plot(angle_spec, 10*log10(x));
title_text = "DOA Angles: " + doa(1) + ", " + doa(2);
title(title_text)

%% CS_Framework_Utilizing_Difference_Coarray 2

CS = CS_Framework_Utilizing_Difference_Coarray(z1, [A1 I2], [1 1]);

CS = CS.Steepest_Descent_DOA();

x = abs(CS.sg(1:end-1));
x = x / max(x);
figure; plot(angle_spec, 10*log10(x));
title_text = "DOA Angles: " + doa(1) + ", " + doa(2);
title(title_text)

%% CS_Framework_Utilizing_Difference_Coarray 3

CS = CS_Framework_Utilizing_Difference_Coarray(z2, A2, [1 1]);

CS = CS.Steepest_Descent_DOA();

x = abs(CS.sg);
x = x / max(x);
figure; plot(angle_spec, 10*log10(x));
title_text = "DOA Angles: " + doa(1) + ", " + doa(2);
title(title_text)