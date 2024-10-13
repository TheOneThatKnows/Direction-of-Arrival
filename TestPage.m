%%

spec = [-inf spec -inf];
[mags, inds] = findpeaks(spec)
spec = spec(2:end-1);

%%
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angles = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
SNR_dB = -10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

spec = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = abs(h' * (y * y') * h);
end
spec = 10 * log10(spec);

plot(angles, spec)
title("CBF")

%%
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angles = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
SNR_dB = -10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

Ry = (1 / L) * (y * y');

spec = DOA.MUSIC(K, 0.5, Ry, sensor_locations, angles);
plot(angles, 10*log10(spec))
title("MUSIC")

%%
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9]; % MRA with 5 sensors
M = length(sensor_locations);
N = sensor_locations(M) + 1;
K = 2;          % # of sources
L = 70;        % # of snapshots

phi_min = 30;
phi_max = 150;
delta_phi = 1;

angles = phi_min:delta_phi:phi_max;

doa = DOA.DOA_Generate(K, phi_min, phi_max, delta_phi);

s = DOA.Source_Generate(K, L);
A = DOA.Array_Manifold(0.5, sensor_locations, doa);
SNR_dB = 10;
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;

Ry = (1 / L) * (y * y');

z = Ry(:);
z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
R_z2 = zeros(N);
for i = 1:N
    z1_i = z1(i:i + N - 1);
    R_z2(:, N-i+1) = z1_i;
end

spec = DOA.MUSIC(K, 0.5, R_z2, 0:N-1, angles);
plot(angles, spec)
title("DA-MUSIC")

%%

t = 30:3:150;
x = sind(t);

figure;
plot(t, x);
xlim([30 150])
grid on

%%
phi_min = 30;
delta_phi = 3;

x = 128.7;

c = floor((x - phi_min) / delta_phi);

phi_a = phi_min + c * delta_phi;
phi_b = phi_a + delta_phi;

percentages = [phi_a phi_b; 1 1] \ [x*100; 100]

%%
clear; clc; close all;

DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);

doa = 60; doa2 = 100; n = 2; % # of sources
a1 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa));
a2 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa2));
N = 100; % # of snapshots
SNR = 0; % dB

s = DOA.Source_Generate(2, N);
v = DOA.Noise_Generate(SNR, M, N);
A = DOA.Array_Manifold(coef, sensor_locations, [doa doa2]);

y = A * s + v;

angles = 0:0.5:180;

% Conventional Beamformer
power_pattern = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    power_pattern(i) = abs(h' * (y * y') * h);
end
power_pattern = power_pattern / max(power_pattern);
% tt = angles;
% [out1, out2] = findpeaks(10 * log10(power_pattern))
figure;
plot(angles, 10 * log10(power_pattern));
hold on

CBF(angles, sensor_locations, y, doa)

% Capon
Ry = y * y'; % covariance matrix
power_pattern2 = zeros(1, length(angles));
for i = 1:length(angles)
    a_ = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    h = (inv(Ry) * a_) / (a_' * inv(Ry) * a_);
    power_pattern2(i) = abs(h' * (y * y') * h);
end
power_pattern2 = power_pattern2 / max(power_pattern2);
plot(angles, 10 * log10(power_pattern2));

% MUSIC
[eig_vecs, ~] = eig(Ry); % eigen decomposition of the covariance matrix
G = eig_vecs(:, 1:M-n);  % noise space
spec = zeros(1, length(angles));
for j = 1:length(angles)
    a_ = exp(1i * pi * sensor_locations.' * cos(angles(j) * pi / 180));
    spec(j) = 1/abs(a_' * (G * G') * a_);
end
spec = spec / max(spec);
plot(angles, 10 * log10(spec))
title("Minimum Redundant Array Processing")
xlim([0 180])
xlabel("phi (degree)")
legend("Conventional", "Capon", "Music")
grid on

%%
clear; clc; close all;
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival');
DOA = FunctionsOfDOA();

s = DOA.Source_Generate(2, 100);
cov(s(1,:), s(2,:))

%%
clear; clc; close all;

L = 100;
n = 0:L-1;
f1 = 200;
f2 = 250;
fs = 4 * max(f1, f2);
% t = 0:6;
x = exp(1i * 2 * pi * (f1 / fs) * n);
y = exp(1i * 2 * pi * (f2 / fs) * n);

v = [var(x); var(y)]
C = cov(x, y)
%% 
% x = zeros(1, length(t));
% x(1:5) = 1;
X = fft(x);
figure;
plot(n, x, '*');
figure;
plot(n, abs(X), '*')

%% 
clear; clc; close all;

w=[0:0.01:2*pi];
k=[0:2*pi/10:2*pi];
Xw = exp(-j*2*w).*sin(5*w/2)./sin(w/2);
Xk = exp(-j*2*k).*sin(5*k/2)./sin(k/2);
plot(w,abs(Xw))
hold
stem(k,abs(Xk),'r')
figure
plot(w,angle(Xw))
hold
stem(k,angle(Xk),'r')
%%
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = [0 1 4 7 9];
M = length(sensor_locations);
N = sensor_locations(M) + 1;

K = 2;
L = 100;
phi_min = 30;
phi_max = 150;
delta_phi = 1;  % angle resolution

doa = (-2 * delta_phi) * ones(1, K);
i = 1;
while true
    temp_angle = phi_min + rand * (phi_max - phi_min);
    temp_array = abs(doa - temp_angle);
    if any(temp_array < delta_phi)
        i = i - 1;
    else
        doa(i) = temp_angle;
    end
    if i == K
        break
    end
    i = i + 1;
end
doa = sort(doa);

doa = round(doa);

A = DOA.Array_Manifold(0.5, sensor_locations, doa);
s = DOA.Source_Generate(K, L);
SNR_dB = 30 * rand;     % min = 0 dB, max = 30 dB
n = DOA.Noise_Generate(SNR_dB, M, L);
y = A * s + n;
R = (1 / L) * (y * y');

z = R(:);
[z1, M_v] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
R_z1 = zeros(M_v);
for i = 1:M_v
    z1_i = z1(i:i + M_v - 1);
    R_z1 = R_z1 + (1 / M_v) * (z1_i * z1_i');
end



%% 
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = 0:63;
M = length(sensor_locations);

doa = 30:0.5:150;
A = DOA.Array_Manifold(0.5, 0:63, doa);
A = DOA.khatri_rao(conj(A), A);
% One = (1 / M) * abs(A(:, 1)' * A(:, 1))
% Zero = (1 / M) * abs(A(:, 1)' * A(:, 2))
I = (1 / M^2) * (A' * A);

L = 1000;
s = DOA.Source_Generate(2, L);
Rs = (1 / L) * (s * s')

%% 
clear; clc; close all;

DOA = FunctionsOfDOA();

sensor_locations = 0:7;
M = length(sensor_locations);

doa = [83 94];
K = length(doa);
SNR_dB = 20;
L = 20; % #of snapshots

A = DOA.Array_Manifold(0.5, sensor_locations, doa);
s = DOA.Source_Generate(K, L);
n = DOA.Noise_Generate(SNR_dB, M, L);

y = A * s + n;
R = (1 / L) * (y * y');

re_R = real(R);
im_R = imag(R);

figure;
imagesc(re_R);
axis equal tight;
colorbar;
title('Re(Original)');

figure;
imagesc(im_R);
axis equal tight;
colorbar;
title('Im(Original)');

%%

sensor_locations = 0:63;
M = length(sensor_locations);

doa = [83 94];
K = length(doa);
SNR_dB = 20;
L = 20; % #of snapshots

A = DOA.Array_Manifold(0.5, sensor_locations, doa);
s = DOA.Source_Generate(K, L);
n = DOA.Noise_Generate(SNR_dB, M, L);

y = A * s + n;
R = (1 / L) * (y * y');

re_R = real(R);
im_R = imag(R);

figure;
imagesc(re_R);
axis equal tight;
colorbar;
title('Re(Desired)');

figure;
imagesc(im_R);
axis equal tight;
colorbar;
title('Im(Desired)');

%%
c = R(:);
R_desired = toeplitz(R);
re_desired = real(R_desired);
im_desired = imag(R_desired);

figure;
imagesc(re_desired);
axis equal tight;
colorbar;
title('Re(Desired)');

figure;
imagesc(im_desired);
axis equal tight;
colorbar;
title('Im(Desired)');
%% 
clear; clc; close all;
load data_6_sensors_3

figure; hold on; grid on;
xlabel('angles (degrees)');
ylabel('Spatial Spectrum (dB)');
title('Spatial Spectra');
plot(angles, 10 * log10(das_ULA));
plot(angles, 10 * log10(capon_ULA));
plot(angles, 10 * log10(music_ULA));
legend('Delay-and-Sum', 'Capon', 'MUSIC');

figure; hold on; grid on;
xlabel('angles (degrees)');
ylabel('Spatial Spectrum (dB)');
title('Spatial Spectra');
plot(angles, 10 * log10(das_ULA));
plot(angles, 10 * log10(capon_ULA));
plot(angles, 10 * log10(music_ULA));
plot(angles, 10 * log10(ss_music_MRA));
legend('Delay-and-Sum', 'Capon', 'MUSIC', 'SS-MUSIC');

%%
clear; clc; close all;

DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

% M = 24;     % # of sensor locations
coef = 0.5; % distance between sensor locations divided by wavelength of signal
% h_loc = zeros(M, 1);
% h_loc([1 2 5 11 17 19 22 24]) = 1; % locations where there are sensors
sensor_locations = [1 2 5 11 17 19 22 24] - 1;
M = length(sensor_locations);

doa = 60; doa2 = 110; n = 2; % # of sources
a1 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa));
a2 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa2));
N = 1000; % # of snapshots
SNR = 10; % dB

s = DOA.Source_Generate(2, N);
v = DOA.Noise_Generate(SNR, M, N);
A = DOA.Array_Manifold(coef, sensor_locations, [doa doa2]);

y = A * s + v;

angles = 0:0.5:180;

% Conventional Beamformer
power_pattern = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    power_pattern(i) = abs(h' * (y * y') * h);
end
power_pattern = power_pattern / max(power_pattern);
figure;
plot(angles, 10 * log10(power_pattern));
hold on

% Capon
Ry = y * y'; % covariance matrix
power_pattern2 = zeros(1, length(angles));
for i = 1:length(angles)
    a_ = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    h = (inv(Ry) * a_) / (a_' * inv(Ry) * a_);
    power_pattern2(i) = abs(h' * (y * y') * h);
end
power_pattern2 = power_pattern2 / max(power_pattern2);
plot(angles, 10 * log10(power_pattern2));

% MUSIC
[eig_vecs, ~] = eig(Ry); % eigen decomposition of the covariance matrix
G = eig_vecs(:, 1:M-n);  % noise space
spec = zeros(1, length(angles));
for j = 1:length(angles)
    a_ = exp(1i * pi * sensor_locations.' * cos(angles(j) * pi / 180));
    spec(j) = 1/abs(a_' * (G * G') * a_);
end
spec = spec / max(spec);
plot(angles, 10 * log10(spec))
title("Minimum Redundant Array Processing")
xlim([0 180])
xlabel("phi (degree)")
legend("Conventional", "Capon", "Music")
grid on

%%
clear; clc; close all;

DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

% M = 24;     % # of sensor locations
coef = 0.5; % distance between sensor locations divided by wavelength of signal
% h_loc = zeros(M, 1);
% h_loc([1 2 5 11 17 19 22 24]) = 1; % locations where there are sensors
sensor_locations = 0:5;
M = length(sensor_locations);

doa = 60; doa2 = 110; n = 2; % # of sources
a1 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa));
a2 = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa2));
N = 1000; % # of snapshots
SNR = 10; % dB

s = DOA.Source_Generate(2, N);
v = DOA.Noise_Generate(SNR, M, N);
A = DOA.Array_Manifold(coef, sensor_locations, [doa doa2]);

y = A * s + v;

angles = 0:0.5:180;

% Conventional Beamformer
power_pattern = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    power_pattern(i) = abs(h' * (y * y') * h);
end
power_pattern = power_pattern / max(power_pattern);
figure;
plot(angles, 10 * log10(power_pattern));
hold on

% Capon
Ry = y * y'; % covariance matrix
power_pattern2 = zeros(1, length(angles));
for i = 1:length(angles)
    a_ = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
    h = (inv(Ry) * a_) / (a_' * inv(Ry) * a_);
    power_pattern2(i) = abs(h' * (y * y') * h);
end
power_pattern2 = power_pattern2 / max(power_pattern2);
plot(angles, 10 * log10(power_pattern2));

% MUSIC
[eig_vecs, ~] = eig(Ry); % eigen decomposition of the covariance matrix
G = eig_vecs(:, 1:M-n);  % noise space
spec = zeros(1, length(angles));
for j = 1:length(angles)
    a_ = exp(1i * pi * sensor_locations.' * cos(angles(j) * pi / 180));
    spec(j) = 1/abs(a_' * (G * G') * a_);
end
spec = spec / max(spec);
plot(angles, 10 * log10(spec))
title("Minimum Redundant Array Processing")
xlim([0 180])
xlabel("phi (degree)")
legend("Conventional", "Capon", "Music")
grid on
%%
clear; clc; close all;

snacd = SNACD(3, 2, 4);
doa = 90 - asind(0:0.2:0.4) - 5;
snacd = snacd.Set_Doa_Angles(doa);
snacd = snacd.Prepare_Array_Manifold();
snacd = snacd.Simulate(10, 1000, true);

% A_p1 = snacd.A_p1;
% Phi_2 = snacd.Phi_2;
% P = snacd.P;

% rhs = A_p1 * Phi_2^(1 - 2*P)
% PI = zeros(2*P);
% for i = 1:2*P
%     PI(i, end-i+1) = 1;
% end
% lhs = PI * conj(A_p1)

% ss = DOA.KR_MUSIC_SNACD(snacd.R, snacd.n, snacd.subarray1_locations, snacd.subarray2_locations, 0.5);
plot(0:0.5:180, 10*log10(snacd.Spatial_Spectrum_Music))

%% CS_Framework
clear; clc; close all;

DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

input = [0 3 3];      % [nested_array M1 M2]
% input = [1 2 3];      % [coprime_array M1 M2]
% input = [2 4 3];      % [super_nested_array M1 M2]
% input = [3 4 4 3];    % [augmented_nested_array_1 M1 M2 L1]
% input = [4 10 10];      % [augmented_nested_array_2 M1 M2]
% input = [5 3 3];      % [nested_array_v2 M1 M2]
% input = [6 2 3 3];    % [sparse_nested_array_with_coprime_displacement_1 N M L]
% input = [7 1 2 3];    % [sparse_nested_array_with_coprime_displacement_2 N M L]

% M: # of sensors
if input(1) == 1 % coprime array
    M = 2 * input(2) - 1 + input(3);
elseif abs(6.5 - input(1)) < 1 % sparse nested array with coprime displacement
    M = 2 * input(4);
else
    M = sum(input(2:3));
end

sensor_locations = DOA.Sensor_Locations(input);             % sensor locations
% sensor_locations = [0 1 4 7 9];
% M = length(sensor_locations);
doa = 90 - asind(0:0.2:0.6) - 5;                            % source angles
% doa = [88 93 165];
snapshots = 100;                                              % # of snapshots
SNR_dB = 15;                                                 % signal to noise ratio in decibels
C = DOA.Mutual_Coupling(200, 0.1, M, sensor_locations);     % mutual coupling

n = length(doa);                                        % number of sources
s = DOA.Source_Generate(n, snapshots);                  % source signals
v = DOA.Noise_Generate(SNR_dB, M, snapshots);           % noise
A = DOA.Array_Manifold(coef, sensor_locations, doa);    % Array Manifold

y = C * A * s + v;
angles = 0:1:180;
A = DOA.Array_Manifold(coef, sensor_locations, angles);

%% On-Grid

CS = CS_Framework(y, A);
CS = CS.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
if snapshots ~= 1
    x = var(CS.Sg.').';
    x = abs(x);
else
    x = abs(CS.Sg);
end
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS1');

Ry = (1/snapshots) * (y * y');
z = Ry(:);
A_d = DOA.khatri_rao(conj(A), A);
I = eye(M); I = I(:);
CS2 = CS_Framework_Utilizing_Difference_Coarray(z, [A_d I], [1 1]);
% CS2 = CS2.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
CS2 = CS2.Steepest_Descent_DOA();
x = abs(CS2.sg(1:end-1));
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS2');

[z1, ~] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
[A1, ~] = DOA.Rearrange_According_to_Sensor_Locations(A_d, sensor_locations);
uDOF = length(z1);
M_v = 0.5 * (uDOF + 1);
I = zeros(uDOF, 1); I(M_v) = 1;
CS3 = CS_Framework_Utilizing_Difference_Coarray(z1, [A1 I], [1 1]);
% CS3 = CS3.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
CS3 = CS3.Steepest_Descent_DOA();
x = abs(CS3.sg(1:end-1));
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS3');

z2 = [z1(1:M_v-1); z1(M_v+1:end)];
A2 = [A1(1:M_v-1, :); A1(M_v+1:end, :)];
CS4 = CS_Framework_Utilizing_Difference_Coarray(z2, A2, [1 2]);
CS4 = CS4.Steepest_Descent_DOA();
x = abs(CS4.sg);
% x = x / max(x);
% figure; plot(angles, 10*log10(x)); title('CS4');
figure; plot(angles, x); title('CS4');

%% Off-Grid

L = -1i * 2 * pi * coef * diag(sind(angles));
D = diag(sensor_locations);
B = D * A * L; % B: Derivative of A

CS_Off = CS_Off_Grid_Framework(y, A, B);
CS_Off = CS_Off.Steepest_Descent_DOA([0.05 0.05]);
if snapshots ~= 1
    x = var(CS_Off.Sg.').';
    x = abs(x);
else
    x = abs(CS_Off.Sg);
end
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS_Off');

%% Gridless-DOA

GD = Gridless_DOA();

T = (1/snapshots) * (y * y');
[GD, doa_angles_ivd, z, c] = GD.IVD(T, sensor_locations, n);
doa_angles_ivd

angles = 0:180;
spec = GD.D2(angles);
figure; plot(angles, abs(spec));
title('IVD using T');

[GD, doa_angles_ap_gridless] = GD.AP_Gridless(y, sensor_locations, n);
doa_angles_ap_gridless

spec = GD.D2(angles);
figure; plot(angles, abs(spec));
title('IVD using Y');

doa

%%

D = 600000;
for i = 1:D
    K = randi(N);
    k = zeros(N, 1);
    k(K) = 1;
    source_angles = (1 + rand(1, K)) * 60;
    source_angles = sort(source_angles);
    for j = 2:K
        if abs(source_angles(j) - source_angles(j-1))
            source_angles(j) = (1 + rand) * 60;
            source_angles = sort(source_angles);
            j = 2;
        end
    end
    A_ohm = Array_Manifold(sensor_locations, source_angles);
    s = Source_Generate(K, L);
    C = Mutual_Coupling(100, 0.1, M, sensor_locations);
    SNR_dB = (rand * 2 - 1) * 15;
    v = Noise_Generate(SNR_dB, M, L);
    y = C * A_ohm * s + v;
    R_ohm = (1/L) * (y * y');
    A = Array_Manifold(0:N-1, source_angles);
    Rs = (1/L) * (s * s');
    R = A * Rs * A';
    u = R(:, 1);

end

D = 600000;
for i = 1:D
    K = randi(N);
    k = zeros(N, 1);
    k(K) = 1;
    source_angles = (1 + rand(1, K)) * 60;
    source_angles = sort(source_angles);
    for j = 2:K
        if abs(source_angles(j) - source_angles(j-1))
            source_angles(j) = (1 + rand) * 60;
            source_angles = sort(source_angles);
            j = 2;
        end
    end
    A_ohm = Array_Manifold(sensor_locations, source_angles);
    s = Source_Generate(K, L);
    C = Mutual_Coupling(100, 0.1, M, sensor_locations);
    SNR_dB = (rand * 2 - 1) * 15;
    v = Noise_Generate(SNR_dB, M, L);
    y = C * A_ohm * s + v;
    R_ohm = (1/L) * (y * y');

    % here we do something that is not important
end










function doa_est = CBF(angles, sensor_locations, y, doa)
spec = zeros(1, length(angles));
for i = 1:length(angles)
    h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
    spec(i) = abs(h' * (y * y') * h);
end
spec = 10 * log10(spec);
[mags, inds] = findpeaks(spec);
doa_est = zeros(2, 1);
[~, ind] = max(mags);
idx = inds(ind);
doa_est(1) = angles(idx);
mags = [mags(1:ind-1) mags(ind+1:end)];
inds = [inds(1:ind-1) inds(ind+1:end)];
[~, ind] = max(mags);
idx = inds(ind);
doa_est(2) = angles(idx);

[~, ind] = min(abs(doa - doa_est(1)));
if ind == 2
    doa = doa(2:-1:1);
end

if rmse(doa, [doa_est(1) doa_est(1)]) < rmse(doa, doa_est)
    doa_est = [doa_est(1) doa_est(1)];
end
end



