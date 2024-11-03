%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Start

% input = [0 7 4];      % [nested_array M1 M2]
% input = [1 2 3];      % [coprime_array M1 M2]
% input = [2 7 4];      % [super_nested_array M1 M2]
% input = [3 4 4 3];    % [augmented_nested_array_1 M1 M2 L1]
% input = [4 7 4];      % [augmented_nested_array_2 M1 M2]
input = [5 3 3];      % [nested_array_v2 M1 M2]
% input = [6 2 3 3];    % [sparse_nested_array_with_coprime_displacement_1 N M L]
% input = [7 1 2 3];    % [sparse_nested_array_with_coprime_displacement_2 N M L]

% M: # of sensors
if input(1) == 1 % coprime array
    M = 2 * input(2) - 1 + input(3)
elseif abs(6.5 - input(1)) < 1 % sparse nested array with coprime displacement
    M = 2 * input(4);
else
    M = sum(input(2:3))
end

%%

% M = 6;
% sensor_locations = 0:5;
% sensor_locations = [0 1 2 6 10 13];

sensor_locations = DOA.Sensor_Locations(input)                  % sensor locations
sensor_placement = DOA.Sensor_Placement(sensor_locations)       % sensor places (in order to visualize the array)
diff_coarray = DOA.Diff_Coarray(sensor_locations)               % difference coarray
uDOF = DOA.Uniform_Degrees_Of_Freedom(sensor_locations)             % uniform degrees of freedom
LuDOF = DOA.One_Side_Uniform_Degrees_Of_Freedom(sensor_locations)   % one side uniform degrees of freedom

%%

% Now, let's use the calculated sensor_locations to compute the array pattern
DOA.Array_Pattern(sensor_locations, coef)

% and solve an example DOA problem
% doa = sort(randi(179, 1, 2));
doa = [50 100 150];
% doa = 90 - asind(0:0.2:0.6) - 5;                            % source angles
snapshots = 1000;                                           % # of snapshots
SNR_dB = 10;                                                 % signal to noise ratio in decibels
C = DOA.Mutual_Coupling(100, 0.1, M, sensor_locations);     % mutual coupling

das_ULA = DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "cbf", C);
capon_ULA = DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "capon", C);
music_ULA = DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "music", C);
DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "kr-music", C);
ss_music_MRA = DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "ss-music", C);
DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "ss-capon", C);