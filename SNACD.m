classdef SNACD < FunctionsOfDOA
    properties
        subarray1_locations;
        subarray2_locations;
        sensor_locations;
        A;      % array manifold matrix
        A1;     % array manifold matrix for subarray 1
        A2;     % array manifold matrix for subarray 2
        N;
        M;
        L;
        doa;    % doa angles of sources
        n;      % # of sources
        y;      % observation matrix
        y1;     % observation matrix for subarray 1
        y2;     % observation matrix for subarray 2
        R;      % cross-covariance matrix
    end
    methods
        function obj = SNACD(N, M, L)
            obj.N = N;
            obj.M = M;
            obj.L = L;
            obj.subarray1_locations = ((1-L)*M:M:0) - N + ((L-1) * M + N);
            obj.subarray2_locations = (0:L*M:(L-1)*L*M) + ((L-1) * M + N);

            obj.sensor_locations = [obj.subarray1_locations obj.subarray2_locations];
        end

        function obj = Set_Doa_Angles(obj, doa)
            obj.doa = doa;
            obj.n = length(doa);
        end

        function obj = Prepare_Array_Manifold(obj)
            obj.A = obj.Array_Manifold(0.5, obj.sensor_locations, obj.doa);
            obj.A1 = obj.A(1:obj.L, :);
            obj.A2 = obj.A(obj.L+1:end, :);
        end

        function obj = Simulate(obj, SNR_dB, N, mutual_coupling)
            s = obj.Source_Generate(obj.n, N);
            v = obj.Noise_Generate(SNR_dB, 2*obj.L, N);

            C = 1;
            if mutual_coupling
                C = obj.Mutual_Coupling(100, 0.1, 2*obj.L, obj.sensor_locations);
            end

            obj.y = C * obj.A * s + v;
            obj.y1 = obj.y(1:obj.L, :);
            obj.y2 = obj.y(obj.L+1:end, :);

            obj.R = (1 / N) * (obj.y1 * obj.y2');
        end
    end
end