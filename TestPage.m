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

% input = [0 10 10];      % [nested_array M1 M2]
% input = [1 2 3];      % [coprime_array M1 M2]
% input = [2 4 3];      % [super_nested_array M1 M2]
% input = [3 4 4 3];    % [augmented_nested_array_1 M1 M2 L1]
input = [4 10 10];      % [augmented_nested_array_2 M1 M2]
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
% doa = 90 - asind(0:0.2:0.6) - 5;                            % source angles
doa = [88 93 165];
snapshots = 1;                                              % # of snapshots
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
CS2 = CS_Framework_Utilizing_Difference_Coarray(z, [A_d I]);
% CS2 = CS2.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
CS2 = CS2.Steepest_Descent_DOA([1 1]);
x = abs(CS2.sg(1:end-1));
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS2');

[z1, ~] = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
[A1, ~] = DOA.Rearrange_According_to_Sensor_Locations(A_d, sensor_locations);
uDOF = length(z1);
M_v = 0.5 * (uDOF + 1);
I = zeros(uDOF, 1); I(M_v) = 1;
CS3 = CS_Framework_Utilizing_Difference_Coarray(z1, [A1 I]);
% CS3 = CS3.Steepest_Descent_DOA([1 min(1/snapshots, 0.1)]);
CS3 = CS3.Steepest_Descent_DOA([1 1]);
x = abs(CS3.sg(1:end-1));
x = x / max(x);
figure; plot(angles, 10*log10(x)); title('CS3');

z2 = [z1(1:M_v-1); z1(M_v+1:end)];
A2 = [A1(1:M_v-1, :); A1(M_v+1:end, :)];
CS4 = CS_Framework_Utilizing_Difference_Coarray(z2, A2);
CS4 = CS4.Steepest_Descent_DOA([1 2]);
x = abs(CS4.sg);
% x = x / max(x);
% figure; plot(angles, 10*log10(x)); title('CS4');
figure; plot(angles, x); title('CS4');

%% Off-Grid

L = -1i * 2 * pi * coef * diag(sind(angles));
D = diag(sensor_locations);
B = D * A * L; % B: Derivative of A

CS_Off = CS_Off_Grid_Framework(y, A, B);
CS_Off = CS_Off.Steepest_Descent_DOA([0.5 0.5]);
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