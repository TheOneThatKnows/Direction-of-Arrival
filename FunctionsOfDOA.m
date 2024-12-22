classdef FunctionsOfDOA
    methods
        function Array_Pattern(~, sensor_locations)
            angles = 0:0.5:180;
            pattern = zeros(1, length(angles));
            for i = 1:length(angles)
                h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                pattern(i) = abs(sum(h));
            end
            pattern = pattern / max(pattern);

            figure;
            plot(angles, pattern);
            xlabel('phi'); xlim([0 180]);
            title('Pattern - Test Array');
            grid on;
        end

        % DOA Generate
        function doa = DOA_Generate(~, number_of_sources, phi_min, phi_max, delta_phi)
            doa = (-2 * delta_phi) * ones(1, number_of_sources);
            i = 1;
            while true
                temp_angle = phi_min + rand * (phi_max - phi_min);
                temp_array = abs(doa - temp_angle);
                if any(temp_array < delta_phi * 2)
                    i = i - 1;
                else
                    doa(i) = temp_angle;
                end
                if i == number_of_sources
                    break
                end
                i = i + 1;
            end
            doa = sort(doa);
        end

        % DOA Estimate
        function doa_est = DOA_Estimate(~, spec, angle_spec, K)
            spec = [-inf spec -inf];
            [peak_mags, peak_inds] = findpeaks(spec);
            doa_est = zeros(1, K);
            [~, sorted_inds] = sort(peak_mags, "descend");
            peak_inds = peak_inds(sorted_inds);
            peak_inds = [peak_inds zeros(1, K-1)];
            idx = peak_inds(1);
            doa_est(1) = angle_spec(idx - 1);

            for i = 2:K
                idx = peak_inds(i);
                if idx == 0
                    doa_est(i:K) = mean(doa_est(1:i-1));
                    break
                else
                    doa_est(i) = angle_spec(idx - 1);
                end
            end

            doa_est = sort(doa_est);
        end

        % Source Generate Old
        function s = Source_Generate_old(~, number_of_sources, number_of_snapshots, vars)
            if nargin == 3
                vars = ones(number_of_sources, 1);
            end
            s = zeros(number_of_sources, number_of_snapshots);
            theta = 180 * rand(number_of_sources, 1);
            for i = 1:number_of_sources
                t_initial = randi(number_of_snapshots);
                t_final = t_initial + number_of_snapshots - 1;
                s(i, :) = sqrt(vars(i)) * exp(1i * pi * (t_initial:t_final) * cosd(theta(i)));
            end
        end

        % Source Generate
        function s = Source_Generate(~, number_of_sources, number_of_snapshots, vars)
            if nargin == 3
                vars = ones(number_of_sources, 1);
            end
            f_max = 9e9;        % 9 GHz
            f_min = 10e3;       % 10 kHz
            fs = 2 * f_max;     % 18 GHz (sampling frequency)
            
            s = zeros(number_of_sources, number_of_snapshots);
            for i = 1:number_of_sources
                f = rand * (f_max - f_min) + f_min;
                t_initial = rand * (fs / f);
                t_final = t_initial + (number_of_snapshots - 1);
                t_range = linspace(t_initial, t_final, number_of_snapshots);
                s(i, :) = sqrt(vars(i)) * exp(1i * 2 * pi * (f / fs) * t_range);
            end
        end

        % Source Generate
        function [s, lambda] = Source_Generate_2(~, number_of_sources, number_of_snapshots, vars)
            if nargin == 3
                vars = ones(number_of_sources, 1);
            end
            f_max = 9e9;        % 9 GHz
            f_min = 10e3;       % 10 kHz
            fs = 2 * f_max;     % 18 GHz (sampling frequency)
            
            s = zeros(number_of_sources, number_of_snapshots);
            lambda = zeros(number_of_sources, 1);
            for i = 1:number_of_sources
                f = rand * (f_max - f_min) + f_min;
                t_initial = rand * (fs / f);
                t_final = t_initial + (number_of_snapshots - 1);
                t_range = linspace(t_initial, t_final, number_of_snapshots);
                s(i, :) = sqrt(vars(i)) * exp(1i * 2 * pi * (f / fs) * t_range);
                lambda(i) = 3e8 / f;
            end
        end

        % Coherent Source Generate
        function s = Coherent_Source_Generate(~, number_of_sources, number_of_snapshots, vars)
            if nargin == 3
                vars = ones(number_of_sources, 1);
            end
            f_max = 9e9;        % 9 GHz
            f_min = 10e3;       % 10 kHz
            fs = 2 * f_max;     % 18 GHz (sampling frequency)
            
            s = zeros(number_of_sources, number_of_snapshots);
            f = rand * (f_max - f_min) + f_min;
            for i = 1:number_of_sources
                t_initial = rand * (fs / f);
                t_final = t_initial + (number_of_snapshots - 1);
                t_range = linspace(t_initial, t_final, number_of_snapshots);
                s(i, :) = sqrt(vars(i)) * exp(1i * 2 * pi * (f / fs) * t_range);
            end
        end

        % Source Generation Using White Noise
        function s = Source_Generate_Ultra(~, number_of_sources, number_of_snapshots, amplitude)
            if nargin == 3
                amplitude = ones(number_of_sources(1), 1);
            end

            f0 = 9e9;        % 9 GHz
            fs = 8 * f0;     % 18 GHz (sampling frequency)
            
            s = zeros(number_of_sources(1), number_of_snapshots);
            for i = 1:number_of_sources(1)
                if i <= number_of_sources(2) && i > 1  % Sources coherent with first source
                    s(i,:) = amplitude(i) * s(1,:);  % Coherent with source 1
                else
                    % Complex Gaussian random process with narrowband filtering
                    white_noise = (randn(1, number_of_snapshots) + 1j*randn(1, number_of_snapshots))/sqrt(2);

                    % Create narrowband filter
                    filter_order = 64;
                    bandwidth = 20; % Hz
                    h = fir1(filter_order, [f0-bandwidth f0+bandwidth]/(fs/2));

                    % Filter the noise to create narrowband signal
                    filtered_signal = filter(h, 1, white_noise);
                    sigma = sqrt(abs(filtered_signal * filtered_signal'));
                    s(i,:) = amplitude(i) * filtered_signal * sqrt(number_of_snapshots) / sigma;
                end
            end
        end

        % Simulate The Environment
        function y = Simulate_Environment(~, sensor_locations, doa, number_of_snapshots, Rs, SNR_dB)

            y = sensorsig(sensor_locations * 0.5, number_of_snapshots, 90-doa, db2pow(-SNR_dB), Rs).';
        end

        % Signal Covariance Matrix
        function Rs = Signal_Covariance(~, K, K_coherent)
            alpha = [1 (0.75 + 0.25 * rand(1, K_coherent - 1))].';
            Rs = eye(K);
            Rs(1:K_coherent, 1:K_coherent) = alpha * alpha';
        end

        % Narrowband Source Generate
        function s = Narrowband_Source_Generate(~, K, L, vars)
            if nargin == 3
                vars = ones(K, 1);
            end
            fc = 9e9;       % 9 GHz
            fs = 2 * fc;    % 18 GHz (sampling frequency)

            v = zeros(K, L);
            for i = 1:K
                v(i, :) = sqrt(vars(i)) * (randn(1, L) + 1i * randn(1, L));
            end
            s = real(v .* exp(1i * (2 * pi * (fc / fs) * (0:L-1))));
        end

        % Narrowband Coherent Source Generate
        function s = Narrowband_Coherent_Source_Generate(~, K, L, vars)
            if nargin == 3
                vars = ones(K, 1);
            end
            fc = 9e9;       % 9 GHz
            fs = 2 * fc;    % 18 GHz (sampling frequency)

            v_og = randn(1, L) + 1i * randn(1, L);
            v = zeros(K, L);
            for i = 1:K
                v(i, :) = sqrt(vars(i)) * v_og;
            end
            s = real(v .* exp(1i * (2 * pi * (fc / fs) * (0:L-1))));
        end

        % Source Generate Final
        function s = Source_Generate_Final(~, K, K_coherent, L, vars, alpha)
            if nargin == 5
                alpha = [1; (0.75+0.25*rand(K_coherent-1, 1)) .* exp(1i * pi * rand(K_coherent-1, 1))];
            end
            v = zeros(K, L);
            v(1:K_coherent, :) = ones(K_coherent, 1) * (randn(1, L) + 1i * randn(1, L));
            v(K_coherent+1:end, :) = randn(K-K_coherent, L) + 1i * randn(K-K_coherent, L);

            s = v .* exp(1i * (2 * pi * 0.5 * (0:L-1)));
            for i = 1:K
                s(i, :) = s(i, :) / sqrt(var(s(i, :)));
            end
            s = [sqrt(vars(1)) * alpha; sqrt(vars(2:end))] .* s;
            shuffledInds = randperm(K);
            s = s(shuffledInds, :);
        end

        % Noise Generate
        function v = Noise_Generate(~, SNR_dB, number_of_sensors, number_of_snapshots)
            v = sqrt(1 / 10^(SNR_dB/10)) * (1/sqrt(2)) * (randn(number_of_sensors, number_of_snapshots) + 1i * randn(number_of_sensors, number_of_snapshots));
        end

        % Array Manifold
        function A = Array_Manifold(~, sensor_locations, doa)
            number_of_sensors = length(sensor_locations);
            number_of_sources = length(doa);
            A = zeros(number_of_sensors, number_of_sources);
            for i = 1:number_of_sources
                A(:, i) = exp(1i * pi * sensor_locations.' * cosd(doa(i)));
            end
        end

        % Array Manifold
        function A = Array_Manifold_2(~, lambda, d, sensor_locations, doa)
            number_of_sensors = length(sensor_locations);
            number_of_sources = length(doa);
            A = zeros(number_of_sensors, number_of_sources);
            for i = 1:number_of_sources
                A(:, i) = exp(1i * 2 * pi * (d / lambda(i)) * sensor_locations.' * cosd(doa(i)));
            end
        end

        % CBF (Delay-and-Sum)
        function spec = CBF(~, y, sensor_locations, angles)
            spec = zeros(1, length(angles));
            for i = 1:length(angles)
                h = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                spec(i) = abs(h' * (y * y') * h);
            end
            spec = spec / max(spec);
        end

        % Capon
        function spec = Capon(~, Ry, sensor_locations, angles)
            spec = zeros(1, length(angles));

            M = size(Ry, 1);

            for i = 1:length(angles)
                a_ = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                h = ((eye(M)/(Ry)) * a_) / (a_' * (eye(M)/(Ry)) * a_);
                spec(i) = abs(h' * Ry * h);
            end
            spec = spec / max(spec);
        end

        % MUSIC
        function spec = MUSIC(~, K, R, sensor_locations, angles)
            spec = zeros(1, length(angles));

            M_v = size(R, 1);

            [eig_vecs, eig_vals] = eig(R, "vector");   % eigen decomposition of the covariance matrix
            [~, sorted_inds] = sort(eig_vals, "ascend");
            eig_vecs = eig_vecs(:, sorted_inds);
            G = eig_vecs(:, 1:M_v-K);   % noise space

            for i = 1:length(angles)
                a_ = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                spec(i) = 1/abs(a_' * (G * G') * a_);
            end
            spec = spec / max(spec);
        end

        % SS-MUSIC
        function [spec, R_z1] = SS_MUSIC(obj, K, R, sensor_locations, angles)
            z = R(:);
            [z1, M_v] = obj.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
            R_z1 = zeros(M_v);
            for i = 1:M_v
                z1_i = z1(i:i + M_v - 1);
                R_z1 = R_z1 + (1 / M_v) * (z1_i * z1_i');
            end
            spec = obj.MUSIC(K, R_z1, 0:M_v-1, angles);
        end

        % DA-MUSIC
        function [spec, R_z2] = DA_MUSIC(obj, K, R, sensor_locations, angles)
            z = R(:);
            [z1, M_v] = obj.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
            R_z2 = zeros(M_v);
            for i = 1:M_v
                z1_i = z1(i:i + M_v - 1);
                R_z2(:, M_v-i+1) = z1_i;
            end
            spec = obj.MUSIC(K, R_z2, 0:M_v-1, angles);
        end

        % DML
        function spec = DML(~, Ry, sensor_locations, angles)
            M = size(Ry, 1);
            spec = zeros(1, length(angles));
            for i = 1:length(angles)
                a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                PI_A_ort = eye(M) - a * (1 / (a' * a)) * a';
                
                spec(i) = trace(PI_A_ort * Ry);
            end
            spec = abs(spec);
            spec = 1 ./ spec;
            spec = spec / max(spec);
        end

        % KR-MUSIC
        function spatial_spectrum = KR_MUSIC(obj, Ry, n, sensor_locations)
            r = Ry(:); % coarray vector
            [S, ~, ~] = svd(r);
            G = S(:, n+1:end); % noise space

            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            for i = 1:length(angles)
                a = exp(1i * pi * sensor_locations.' * cosd(angles(i)));
                kr_product = obj.khatri_rao(conj(a), a);
                spatial_spectrum(i) = 1/abs(kr_product' * (G * G') * kr_product);
            end
            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);
        end
        
        % KR-MUSIC for Sparse Nested Arrays with Coprime Displacement
        function spatial_spectrum = KR_MUSIC_SNACD(obj, Ry, n, subarray1_locations, subarray2_locations)
            r = Ry(:); % coarray vector
            [S, ~, ~] = svd(r);
            G = S(:, n+1:end); % noise space

            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            for i = 1:length(angles)
                a1 = exp(1i * pi * subarray1_locations.' * cosd(angles(i)));
                a2 = exp(1i * pi * subarray2_locations.' * cosd(angles(i)));
                kr_product = obj.khatri_rao(conj(a2), a1);
                spatial_spectrum(i) = 1/abs(kr_product' * (G * G') * kr_product);
            end
            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);
        end

        % Sensor Locations
        function sensor_locations = Sensor_Locations(obj, input)
            if input(1) == 0 % nested array
                sensor_locations = obj.Nested_Array_Locations(input(2:end));

            elseif input(1) == 1 % coprime array
                sensor_locations = obj.Coprime_Array_Locations(input(2), input(3));

            elseif input(1) == 2 % super-nested array
                sensor_locations = obj.Super_Nested_Array_Locations(input(2), input(3));

            elseif input(1) == 3 % augmented-nested array 1
                sensor_locations = obj.Augmented_Nested_Array_1_Locations(input(2), input(3), input(4));

            elseif input(1) == 4 % augmented-nested array 2
                sensor_locations = obj.Augmented_Nested_Array_2_Locations(input(2), input(3));

            elseif input(1) == 5 % nested-array-v2
                sensor_locations = obj.Nested_Array_v2_Locations(input(2), input(3));

            elseif input(1) == 6 % sparse nested array with coprime displacement 1
                sensor_locations = obj.SNACD_1_Locations(input(2), input(3), input(4));

            elseif input(1) == 7 % sparse nested array with coprime displacement 2
                sensor_locations = obj.SNACD_2_Locations(input(2), input(3), input(4));
            end
        end
        
        % Nested Array Locations
        function sensor_locations = Nested_Array_Locations(~, N)
            K = length(N); % number of levels

            sensor_locations = zeros(1, sum(N));
            sensor_locations(1:N(1)) = 1:N(1);
            for i = 2:K
                locations = 1:N(i);
                for j = 1:i-1
                    for n = 1:N(i)
                        locations(n) = locations(n) * (N(j) + 1);
                    end
                end
                idx_init = sum(N(1:i-1))+1;
                idx_end = idx_init + N(i) - 1;
                sensor_locations(idx_init:idx_end) = locations;
                sensor_locations = sensor_locations - 1;
            end
        end
        
        % Nested Array v2 Locations
        function sensor_locations = Nested_Array_v2_Locations(~, N1, N2)
            sensor_locations = zeros(1, N1 + N2);
            sensor_locations(1:N1+1) = 0:N1;
            sensor_locations(N1+2:N1+N2) = (2*(N1 + 1)):(N1 + 1):(N2*(N1+1));
        end
        
        % Coprime Array Locations
        function sensor_locations = Coprime_Array_Locations(~, M, N)
            sensor_locations_1 = M * (0:N-1);
            sensor_locations_2 = N * (0:2*M-1);

            sensor_locations = sort([sensor_locations_1 sensor_locations_2(2:end)]);
        end
        
        % Super-Nested Array Locations
        function sensor_locations = Super_Nested_Array_Locations(~, N1, N2)
            REM = rem(N1, 4);
            r = floor(N1 / 4);
            if REM == 0
                A1 = r; B1 = r - 1; A2 = r - 1; B2 = r - 2;
            elseif REM == 1
                A1 = r; B1 = r - 1; A2 = r - 1; B2 = r - 1;
            elseif REM == 2
                A1 = r + 1; B1 = r - 1; A2 = r; B2 = r - 2;
            elseif REM == 3
                A1 = r; B1 = r; A2 = r; B2 = r - 1;
            end
            X1 = 1 + 2 * (0:A1);
            Y1 = N1 + 1 - (1 + 2 * (0:B1));
            X2 = (N1 + 1) + (2 + 2 * (0:A2));
            Y2 = 2 * (N1 + 1) - (2 + 2 * (0:B2));
            Z1 = (2:N2) * (N1 + 1);
            Z2 = N2 * (N1 + 1) - 1;

            sensor_locations = sort([X1 Y1 X2 Y2 Z1 Z2]);
            sensor_locations = sensor_locations - 1;
        end
        
        % Augmented Nested Array 1 Positions
        function sensor_locations = Augmented_Nested_Array_1_Locations(obj, N1, N2, L1)
            if L1 == -1
                L1 = ceil(0.5 * (N1 + 1));
            end
            temp_locations = obj.Nested_Array_Locations([N1 N2]) + 1; % parent nested array
            temp_locations(N1 - L1 + 1:N1) = (N1 + 1) * N2 + temp_locations(N1 - L1 + 1:N1) - temp_locations(N1 - L1 + 1) + 1;
            sensor_locations = sort(temp_locations) - 1;
        end
        
        % Augmented Nested Array 2 Positions
        function sensor_locations = Augmented_Nested_Array_2_Locations(obj, N1, N2)
            temp_locations = obj.Nested_Array_Locations([N1 N2]) + 1; % parent nested array
            temp_locations(3:2:N1) = (N1 + 1) * N2 - (rem(N1, 2) + 1) + temp_locations(3:2:N1);
            sensor_locations = sort(temp_locations) - 1;
        end

        % Sparse Nested Array with Coprime Displacement 1 Positions
        function sensor_locations = SNACD_1_Locations(~, N, M, L)
            subarray1 = (0:M:(L-1)*M) + N;
            subarray2 = 0:L*M:(L-1)*L*M;

            sensor_locations = sort([subarray1 subarray2]);
        end

        % Sparse Nested Array with Coprime Displacement 2 Positions
        function sensor_locations = SNACD_2_Locations(~, N, M, L)
            subarray1 = ((1-L)*M:M:0) - N;
            subarray2 = 0:L*M:(L-1)*L*M;

            sensor_locations = [subarray1 subarray2] + ((L-1) * M + N);
        end
        
        % Sensor Placement
        function sensor_placement = Sensor_Placement(~, sensor_locations)
            sensor_placement = zeros(1, sensor_locations(end)+1);
            sensor_placement(sensor_locations + 1) = 1;
        end
        
        % Difference Coarray
        function diff_coarray = Diff_Coarray(obj, sensor_locations)
            sensor_placement = obj.Sensor_Placement(sensor_locations);
            N = length(sensor_placement);
            diff_coarray = zeros(1, 2 * N - 1);
            for i = 1:N
                diff_coarray(N + i - 1) = sensor_placement(i:end) * sensor_placement(1:end-i+1).';
            end
            diff_coarray(1:N-1) = diff_coarray(end:-1:N+1);
        end
        
        % Degrees of Freedom
        function degrees_of_freedom = Degrees_Of_Freedom(~, N)
            K = length(N); % number of levels
            degrees_of_freedom = 2 * (([0 N] * [N 0].') + N(K) - 1) + 1;
        end
        
        % Uniform Degrees of Freedom
        function uniform_degrees_of_freedom = Uniform_Degrees_Of_Freedom(obj, sensor_locations)
            diff_coarray = obj.Diff_Coarray(sensor_locations);
            uniform_degrees_of_freedom = 0;
            length_diff_coarray = length(diff_coarray);
            for i = (length_diff_coarray + 1)/2:length_diff_coarray
                if diff_coarray(i) == 0
                    break
                end
                uniform_degrees_of_freedom = uniform_degrees_of_freedom + 1;
            end
            uniform_degrees_of_freedom = 2 * uniform_degrees_of_freedom - 1;
        end
        
        % One-Side Uniform Degrees of Freedom
        function one_side_uniform_degrees_of_freedom = One_Side_Uniform_Degrees_Of_Freedom(obj, sensor_locations)
            uniform_degrees_of_freedom = obj.Uniform_Degrees_Of_Freedom(sensor_locations);
            one_side_uniform_degrees_of_freedom = 0.5 * (uniform_degrees_of_freedom - 1);
        end

        % Mutual Coupling
        function C = Mutual_Coupling(~, B, c_s, M, sensor_locations)
            C = eye(M);

            if B < 1 % if mutual_coupling does not exist
                return
            end

            c = zeros(B-1, 1);
            c(1) = c_s * exp(1i * rand * 2 * pi);
            c(2:B-1) = c(1) * exp(-1i * pi * (1:B-2).' / 8) ./ (2:B-1).';
            for i = 1:M-1
                for j = i+1:M
                    lag = abs(sensor_locations(i) - sensor_locations(j));
                    if lag <= B
                        C(i, j) = c(lag);
                        C(j, i) = c(lag);
                    end
                end
            end
        end

        % Functions Used In Spatial Smoothing
        % Khatri-Rao Product
        function X = khatri_rao(~, A1, A2)
            [r1, c] = size(A1);
            [r2, ~] = size(A2);
            X = zeros(r1*r2, c);
            for j = 1:c
                for i = 1:r1
                    row_init = (i-1) * r2 + 1;
                    row_end = i * r2;
                    X(row_init:row_end, j) = A1(i, j) * A2(:, j);
                end
            end
        end

        % Kronecker Product
        function X = kronecker(~, A1, A2)
            [r1, c1] = size(A1);
            [r2, c2] = size(A2);
            X = zeros(r1*r2, c1*c2);
            for i = 1:r1
                for j = 1:c1
                    row_interval = (i-1)*r2+1:i*r2;
                    col_interval = (j-1)*c2+1:j*c2;
                    X(row_interval, col_interval) = A1(i, j) * A2;
                end
            end
        end

        % Sort and Discard Repeating Rows According to Sensor Locations In The Difference Coarray
        function [X2, M_v] = Rearrange_According_to_Sensor_Locations(obj, X1, sensor_locations)
            diff_coarray = obj.Diff_Coarray(sensor_locations);
            uDOF = obj.Uniform_Degrees_Of_Freedom(sensor_locations);
            M_v = obj.One_Side_Uniform_Degrees_Of_Freedom(sensor_locations) + 1;
            M = length(sensor_locations);

            [~, c] = size(X1);
            X2 = zeros(uDOF, c);
            diff_vector = zeros(uDOF, 1);
            for i = 1:length(sensor_locations)
                for j = 1:length(sensor_locations)
                    diff = -sensor_locations(i) + sensor_locations(j);
                    idx2 = M_v + diff;
                    if (idx2 < 1 || idx2 > uDOF) || diff_vector(idx2) == diff_coarray(idx2)
                        continue
                    end

                    idx1 = (i-1) * M + j;
                    X2(idx2, :) = X2(idx2, :) + (1/diff_coarray(idx2)) * X1(idx1, :);
                    diff_vector(idx2) = diff_vector(idx2) + 1;
                end
            end
        end
    end
end