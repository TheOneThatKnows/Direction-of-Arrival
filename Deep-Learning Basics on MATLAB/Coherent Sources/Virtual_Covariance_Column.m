function column = Virtual_Covariance_Column(DOA, R, sensor_locations)
M = length(sensor_locations);
N = sensor_locations(M) - sensor_locations(1) + 1;
diff_coarray = DOA.Diff_Coarray(sensor_locations);

sparse_entries = zeros(2*N-1, 1);
for i = 1:M
    for j = 1:M
        diff = sensor_locations(i) - sensor_locations(j);
        count = diff_coarray(N + diff);
        sparse_entries(N + diff) = sparse_entries(N + diff) + (1 / count) * R(i, j);
    end
end
column = [sparse_entries(N); 0.5 * (sparse_entries(N+1:1:2*N-1) + conj(sparse_entries(N-1:-1:1)))];
end
