%% Initialization

clear; clc; close all;
DOA = FunctionsOfDOA();
coef = 0.5; % unit distance between sensors divided by the wavelength of the signal

%% Start

input = [0 3 3];        % [nested_array M1 M2]
% input = [1 2 3];      % [coprime_array M1 M2]
% input = [2 4 3];      % [super_nested_array M1 M2]
% input = [3 4 4 3];    % [augmented_nested_array_1 M1 M2 L1]
% input = [4 7 4];      % [augmented_nested_array_2 M1 M2]
% input = [5 3 3];        % [nested_array_v2 M1 M2]

if input(1) ~= 1
    M = sum(input(2:3))                                     % # of sensors
else
    M = 2 * input(2) - 1 + input(3)
end
sensor_locations = DOA.Sensor_Locations(input)                  % sensor locations
sensor_placement = DOA.Sensor_Placement(sensor_locations)       % sensor places (in order to visualize the array)
diff_coarray = DOA.Diff_Coarray(sensor_placement)               % difference coarray
uDOF = DOA.Uniform_Degrees_Of_Freedom(diff_coarray)             % uniform degrees of freedom
LuDOF = DOA.One_Side_Uniform_Degrees_Of_Freedom(diff_coarray)   % one side uniform degrees of freedom

% Now, let's use the calculated sensor_locations to compute the array pattern
DOA.Array_Pattern(sensor_locations, coef)

% and solve an example DOA problem
% doa = sort(randi(179, 1, 2));
% doa = [50 150];
doa = 90 - asind(0:0.2:0.6) - 5;                            % source angles
snapshots = 1000;                                           % # of snapshots
SNR_dB = 0;                                                 % signal to noise ratio in decibels
C = DOA.Mutual_Coupling(100, 0.1, M, sensor_locations);     % mutual coupling

DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "cbf", C);
DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "capon", C);
DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "music", C);
DOA.Spatial_Spectrum(M, sensor_locations, doa, snapshots, SNR_dB, coef, "kr-music", C);