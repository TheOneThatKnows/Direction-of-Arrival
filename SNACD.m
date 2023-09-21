classdef SNACD < FunctionsOfDOA
    properties
        subarray1_locations;
        subarray2_locations;
        sensor_locations;
        diff_vector;
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
        r;      % vectorized version of R
        P;      % max integer smaller than (L^2)/3
        Q;      % number of snapshots
        Y;      % matrix whose columns are q-th to (q + 2P-1)-th row of r
        A_p1;   % direction matrix of r_1
        Phi_2;  % diag(exp(-1i *  pi * M * cosd(doa)))
        B;      % Y = A_p1 * B
        PI_2P;
        PI_Q;
        Z;      % augmented output
        Bz;     % Z = A_p1 * Bz
        Rz;     % covariance matrix
        Spatial_Spectrum_Music;
    end
    methods
        function obj = SNACD(N, M, L)
            obj.N = N;
            obj.M = M;
            obj.L = L;
            
            obj.P = floor((L^2)/3) - 0.5 * (1 - sign(rem(L, 3) - 0.5));
            obj.Q = L^2 - 2 * obj.P + 1;
            
            obj.subarray1_locations = ((1-L)*M:M:0);
            obj.subarray2_locations = (0:L*M:(L-1)*L*M);

            obj.sensor_locations = [obj.subarray1_locations obj.subarray2_locations];

            obj.diff_vector = obj.Diff_Vector(obj.subarray1_locations, obj.subarray2_locations);

            obj.PI_2P = obj.PI(2 * obj.P);
            obj.PI_Q = obj.PI(obj.Q);
        end

        function obj = Set_Doa_Angles(obj, doa)
            obj.doa = doa;
            obj.n = length(doa);
            obj.Phi_1 = diag(exp(-1i * pi * obj.N * cosd(doa)));
            obj.Phi_2 = diag(exp(-1i * pi * obj.M * cosd(doa)));
        end

        function obj = Prepare_Array_Manifold(obj)
            obj.A1 = obj.Array_Manifold(0.5, obj.subarray1_locations, obj.doa);
            obj.A2 = obj.Array_Manifold(0.5, obj.subarray2_locations, obj.doa);
            A_pre = obj.khatri_rao(conj(obj.A2), obj.A1);
            obj.A = obj.Order_By_Diff_Vector(A_pre, obj.diff_vector);
            obj.A_p1 = obj.A(1:2*obj.P, :);
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
            r_pre = obj.R(:);
            obj.r = obj.Order_By_Diff_Vector(r_pre, obj.diff_vector);
            
            vars = var(s.').';
            obj.mu = obj.Phi_1 * vars;

            obj.Y = zeros(2*obj.P, obj.Q);
            for q = 1:obj.Q
                obj.Y(:, q) = obj.r(q:(q + 2 * obj.P - 1));
            end
            obj.B = zeros(obj.n, obj.Q);
            for i = 1:obj.Q
                obj.B(:, i) = obj.Phi_2^(i-1) * obj.mu;
            end

            obj.Z = [obj.Y (obj.PI_2P * conj(obj.Y) * obj.PI_Q)];
            obj.Bz = [obj.B (obj.Phi_2^(1 - 2*obj.P) * conj(obj.B) * obj.PI_Q)];

            obj.Rz = (1 / (2*obj.L)) * (obj.Z * obj.Z');

            obj.Spatial_Spectrum_Music = obj.MUSIC_SNACD();
        end

        function v = Diff_Vector(obj, v1, v2)
            v = zeros(1, obj.L^2);
            for i = 1:obj.L
                start_idx = (i - 1) * obj.L + 1;
                end_idx = start_idx + obj.L - 1;

                v(start_idx:end_idx) = -v2(i) + v1;
            end
        end

        function X_out = Order_By_Diff_Vector(obj, X_in, diff_vector)
            [row, column] = size(X_in);
            X_out = zeros(row, column);
            for i = 1:obj.L^2
                diff_val = 0 - (i - 1) * obj.M;
                [~, idx] = min(abs(diff_vector - diff_val));
                X_out(i, :) = X_in(idx, :);
            end
        end

        function Pi = PI(~, num)
            Pi = zeros(num);
            for i = 1:num
                Pi(i, end-i+1) = 1;
            end
        end

        function spatial_spectrum = MUSIC_SNACD(obj)
            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            [eig_vecs, ~] = eig(obj.Rz);   % eigen decomposition of the covariance matrix
            G = eig_vecs(:, 1:(2*obj.P-obj.n));   % noise space

            for i = 1:length(angles)
                a_p1 = exp(1i * pi * (0:(2*obj.P-1)).' * obj.M * cosd(angles(i)));
                spatial_spectrum(i) = 1/abs(a_p1' * (G * G') * a_p1);
            end
            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);
        end
    end
end