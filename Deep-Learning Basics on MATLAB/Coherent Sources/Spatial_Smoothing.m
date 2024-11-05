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