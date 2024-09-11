classdef Network_Functions
    methods
        function R = CRN2_Function_v1(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 2);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);

            R_c1 = predict(net, feature).';

            N = (length(R_c1) + 1) * 0.5;
            R_c1 = R_c1(1:N) + 1i * [0; R_c1(N+1:end)];
            R = toeplitz(R_c1');
        end

        function spec = CRN2_Function_v2_1(~, crn_network_v2, crn_network_v2_1, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 2);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);

            R_c1 = predict(crn_network_v2, feature(:, :, 1), feature(:, :, 2)).';
            spec = predict(crn_network_v2_1, R_c1.');
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_2(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 2);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);

            spec = predict(net, feature(:, :, 1), feature(:, :, 2));
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_4(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 2);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);

            spec = predict(net, feature(:, :, 1), feature(:, :, 2), feature);
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_5(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 3);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);
            feature(:, :, 3) = angle(R_ohm);

            spec = predict(net, feature(:, :, 1), feature(:, :, 2), feature(:, :, 3));
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_6(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 3);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);
            feature(:, :, 3) = angle(R_ohm) / pi;

            spec = predict(net, feature);
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_7(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 2);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);

            spec = predict(net, feature(:, :, 1), feature(:, :, 2));
            spec = spec / max(spec);
        end

        function spec = CRN2_Function_v2_8(~, net, M, R_ohm)
            normalized_R_ohm = R_ohm / max(diag(abs(R_ohm)));

            feature = zeros(M, M, 3);
            feature(:, :, 1) = real(normalized_R_ohm);
            feature(:, :, 2) = imag(normalized_R_ohm);
            feature(:, :, 3) = angle(R_ohm) / pi;

            spec = predict(net, feature);
            spec = spec / max(spec);
        end

        function spec = SparseDNN_Function_v2(~, net, M, R_ohm)
            feature = zeros(M*M, 1);
            re_R = real(R_ohm);
            im_R = imag(R_ohm);

            feature(1:M) = diag(re_R);
            ind = M + 1;
            for i = 2:M
                for j = 1:i-1
                    feature(ind:ind+1) = [re_R(i, j); im_R(i, j)];
                    ind = ind + 2;
                end
            end

            spec = predict(net, feature.');
            spec = spec / max(spec);
        end

        function spec = SparseDNN_Function_v3(~, net, N, R_ohm, DOA, sensor_locations)
            z = R_ohm(:);
            z1 = DOA.Rearrange_According_to_Sensor_Locations(z, sensor_locations);
            R_z1 = zeros(N);
            for i = 1:N
                z1_i = z1(i:i + N - 1);
                R_z1 = R_z1 + (1 / N) * (z1_i * z1_i');
            end

            re_R_z1 = real(R_z1(:, 1));
            im_R_z1 = imag(R_z1(:, 1));

            feature = [re_R_z1; im_R_z1(2:end)];
            feature = feature / max(feature);

            spec = predict(net, feature.');
            spec = spec / max(spec);
        end

        function spec = Sparse_1D_Function(~, net, Ry, DOA, A_sparse)
            z = Ry(:);
            A1 = DOA.khatri_rao(conj(A_sparse), A_sparse);

            feature = abs(A1' * z);
            feature = (feature - mean(feature)) / var(feature);
            spec = predict(net, feature);
        end
    end
end