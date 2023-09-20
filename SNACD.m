classdef SNACD < FunctionsOfDOA
    properties
        subarray1_locations;
        subarray2_locations;
        sensor_locations;
        A;      % virtual direction matrix
        A1;     % array manifold matrix for subarray 1
        Phi_1;  % translation matrix
        A2;     % array manifold matrix for subarray 2
        mu;     % source vector
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
            obj.subarray1_locations = ((1-L)*M:M:0);
            obj.subarray2_locations = (0:L*M:(L-1)*L*M);

            obj.sensor_locations = [obj.subarray1_locations obj.subarray2_locations];
        end

        function obj = Set_Doa_Angles(obj, doa)
            obj.doa = doa;
            obj.n = length(doa);
            obj.Phi_1 = diag(exp(-1i * pi * obj.N * cosd(doa)));
        end

        function obj = Prepare_Array_Manifold(obj)
            obj.A1 = obj.Array_Manifold(0.5, obj.subarray1_locations, obj.doa);
            obj.A2 = obj.Array_Manifold(0.5, obj.subarray2_locations, obj.doa);
            obj.A = obj.khatri_rao(conj(obj.A2), obj.A1);
        end

        function obj = Simulate(obj, SNR_dB, N, mutual_coupling, vars)
            if nargin == 4
                vars = ones(obj.n, 1);
            end

            s = obj.Source_Generate(obj.n, N, vars);
            v1 = obj.Noise_Generate(SNR_dB, obj.L, N);
            v2 = obj.Noise_Generate(SNR_dB, obj.L, N);

            C1 = 1; C2 = 1;
            if mutual_coupling
                C1 = obj.Mutual_Coupling(100, 0.1, obj.L, obj.subarray1_locations);
                C2 = obj.Mutual_Coupling(100, 0.1, obj.L, obj.subarray2_locations);
            end

            obj.y1 = C1 * obj.A1 * obj.Phi_1 * s + v1;
            obj.y2 = C2 * obj.A2 * s + v2;

            obj.R = (1 / N) * (obj.y1 * obj.y2');
            vars = var(s.').';
            obj.mu = obj.Phi_1 * vars;
        end
    end
end