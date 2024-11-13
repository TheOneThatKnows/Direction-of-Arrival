function R_out = Sparse_Toeplitz(DOA, R, sensor_locations)
M = length(sensor_locations);
N = sensor_locations(M) - sensor_locations(1) + 1;
diff_coarray = DOA.Diff_Coarray(sensor_locations);

sparse_entries = zeros(1, 2*N-1);
for i = 1:M
    for j = 1:M
        diff = sensor_locations(i) - sensor_locations(j);
        count = diff_coarray(N + diff);
        sparse_entries(N + diff) = sparse_entries(N + diff) + (1 / count) * R(i, j);
    end
end

R_out = zeros(M);
for i = 1:M
    for j = 1:M
        diff = sensor_locations(i) - sensor_locations(j);
        R_out(i, j) = sparse_entries(N + diff);
    end
end
end