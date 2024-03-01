classdef FunctionsOfDOA
    methods
        function Array_Pattern(~, sensor_locations, coef)
            angles = 0:0.5:180;
            pattern = zeros(1, length(angles));
            for i = 1:length(angles)
                h = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
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

        % Source Generate
        function s = Source_Generate(~, number_of_sources, number_of_snapshots, vars)
            if nargin == 3
                vars = ones(number_of_sources, 1);
            end
            s = zeros(number_of_sources, number_of_snapshots);
            theta = 180 * rand(number_of_sources, 1);
            for i = 1:number_of_sources
                s(i, :) = sqrt(vars(i)) * exp(1i * pi * (0:number_of_snapshots-1) * cosd(theta(i)));
            end
        end

        % Noise Generate
        function v = Noise_Generate(~, SNR_dB, number_of_sensors, number_of_snapshots)
            v = (1 / 10^(SNR_dB/10)) * (1/sqrt(2)) * (randn(number_of_sensors, number_of_snapshots) + 1i * randn(number_of_sensors, number_of_snapshots));
        end

        % Array Manifold
        function A = Array_Manifold(~, coef, sensor_locations, doa)
            number_of_sensors = length(sensor_locations);
            number_of_sources = length(doa);
            A = zeros(number_of_sensors, number_of_sources);
            for i = 1:number_of_sources
                A(:, i) = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(doa(i)));
            end
        end

        % Spatial Spectrum Computation
        function spatial_spectrum = Spatial_Spectrum(obj, M, sensor_locations, doa, N, SNR, coef, method, C)
            n = length(doa);

            s = obj.Source_Generate(n, N);
            v = obj.Noise_Generate(SNR, M, N);
            A = obj.Array_Manifold(coef, sensor_locations, doa);

            y = C * A * s + v;

            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            if method == "cbf"
                for i = 1:length(angles)
                    h = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
                    spatial_spectrum(i) = abs(h' * (y * y') * h);
                end
            else
                Ry = (1 / N) * (y * y');

                if method == "capon"
                    spatial_spectrum = obj.Capon(coef, Ry, sensor_locations);

                elseif method == "music"
                    spatial_spectrum = obj.MUSIC(n, coef, Ry, sensor_locations);

                elseif method == "kr-music"
                    spatial_spectrum = obj.KR_MUSIC(Ry, n, sensor_locations, coef);

                else % spatial smoothing methods
                    z = Ry(:);
                    [z1, M_v] = obj.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
                    R_z1 = zeros(M_v);
                    for i = 1:M_v
                        z1_i = z1(i:i + M_v - 1);
                        R_z1 = R_z1 + (1 / M_v) * (z1_i * z1_i');
                    end

                    if method == "ss-music"
                        spatial_spectrum = obj.MUSIC(n, coef, R_z1);
                    else % if method == "ss-capon"
                        spatial_spectrum = obj.Capon(coef, R_z1);
                    end
                end
            end

            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);

            figure;
            plot(angles, 10 * log10(spatial_spectrum));
            xlabel('phi'); xlim([0 180]);
            title_str = 'Spatial Spectrum - ' + method;
            title(title_str);
            grid on;
        end

        % Capon
        function spatial_spectrum = Capon(~, coef, Ry, sensor_locations)
            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            M = length(Ry(:, 1));

            if nargin == 3
                sensor_locations = 0:M-1;
            end

            for i = 1:length(angles)
                a_ = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
                h = ((eye(M)/(Ry)) * a_) / (a_' * (eye(M)/(Ry)) * a_);
                spatial_spectrum(i) = abs(h' * Ry * h);
            end
        end

        % MUSIC
        % This function works for a given covariance matrix. I try to use it when I calculate a covariance matrix using spatial smoothing.
        function [spatial_spectrum, angles] = MUSIC(~, n, coef, Rz1, sensor_locations)
            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            M_v = length(Rz1(:, 1));

            [eig_vecs, ~] = eig(Rz1);   % eigen decomposition of the covariance matrix
            G = eig_vecs(:, 1:M_v-n);   % noise space

            if nargin == 4
                sensor_locations = 0:M_v-1; % uniform linear array
            end

            for i = 1:length(angles)
                a_ = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
                spatial_spectrum(i) = 1/abs(a_' * (G * G') * a_);
            end
            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);
        end

        % KR-MUSIC
        function spatial_spectrum = KR_MUSIC(obj, Ry, n, sensor_locations, coef)
            r = Ry(:); % coarray vector
            [S, ~, ~] = svd(r);
            G = S(:, n+1:end); % noise space

            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            for i = 1:length(angles)
                a = exp(1i * 2 * pi * coef * sensor_locations.' * cosd(angles(i)));
                kr_product = obj.khatri_rao(conj(a), a);
                spatial_spectrum(i) = 1/abs(kr_product' * (G * G') * kr_product);
            end
            spatial_spectrum = spatial_spectrum / max(spatial_spectrum);
        end
        
        % KR-MUSIC for Sparse Nested Arrays with Coprime Displacement
        function spatial_spectrum = KR_MUSIC_SNACD(obj, Ry, n, subarray1_locations, subarray2_locations, coef)
            r = Ry(:); % coarray vector
            [S, ~, ~] = svd(r);
            G = S(:, n+1:end); % noise space

            angles = 0:0.5:180;
            spatial_spectrum = zeros(1, length(angles));

            for i = 1:length(angles)
                a1 = exp(1i * 2 * pi * coef * subarray1_locations.' * cosd(angles(i)));
                a2 = exp(1i * 2 * pi * coef * subarray2_locations.' * cosd(angles(i)));
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
        function diff_coarray = Diff_Coarray(~, sensor_placement)
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
        function uniform_degrees_of_freedom = Uniform_Degrees_Of_Freedom(~, diff_coarray)
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
        function one_side_uniform_degrees_of_freedom = One_Side_Uniform_Degrees_Of_Freedom(obj, diff_coarray)
            uniform_degrees_of_freedom = obj.Uniform_Degrees_Of_Freedom(diff_coarray);
            one_side_uniform_degrees_of_freedom = 0.5 * (uniform_degrees_of_freedom - 1);
        end

        % Mutual Coupling
        function C = Mutual_Coupling(~, B, c_s, M, sensor_locations)
            C = eye(M);

            if B < 1 % if mutual_coupling does not exist
                return
            end

            c = zeros(B-1, 1);
            c(1) = c_s * exp(1i * pi / 3);
            c(2:B-1) = c(1) * exp(-1i * pi * (1:B-2).' / 8);
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

        % Sort and Discard Repeating Rows According to Sensor Locations In The Difference Coarray
        function [X2, M_v] = Rearrange_According_to_Sensor_Locations(obj, X1, sensor_locations)
            sensor_placement = obj.Sensor_Placement(sensor_locations);
            diff_coarray = obj.Diff_Coarray(sensor_placement);
            uDOF = obj.Uniform_Degrees_Of_Freedom(diff_coarray);
            M_v = obj.One_Side_Uniform_Degrees_Of_Freedom(diff_coarray) + 1;
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