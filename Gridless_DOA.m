classdef Gridless_DOA
    properties
        G;
        gamma;
        doa_angles;
        z;
        c;
        spectrum;
    end
    methods
        function obj = Gridless_DOA()
        end

        % Irregular Vandermonde Decomposition
        function obj = IVD(obj, T, gamma, K)
            % T: Autocovarience Matrix
            % gamma: Sensor Positions
            % K: # of sources

            obj.gamma = gamma;
            M = length(gamma);  % # of sensors
            [U, ~, ~] = svd(T);
            U_N = U(:, K+1:M);  % Noise Space

            [obj, doa_angles] = obj.estimate_doas(U_N, M, K);
            z = exp(1i * pi * cosd(doa_angles));
            W = zeros(M, K);
            for i = 1:K
                W(:, i) = z(i).^(gamma.');
            end
            W_pseudo_inv = (W' * W)\W';
            C = W_pseudo_inv * T * W_pseudo_inv';
            c = zeros(K, 1);
            for i = 1:K
                c(i) = C(i, i);
            end
            obj.doa_angles = doa_angles;
            obj.z = z;
            obj.c = c;
        end

        function [obj, doa_angles] = estimate_doas(obj, U_N, M, K)
            obj.G = U_N * U_N';
            intervals = 18/M * (0:10*M-1).' * [1 1] + [0 18/M];
            idx = 1;
            doa_candidates = [];
            mags = [];
            for i = 1:10*M
                a = intervals(i, 1);
                b = intervals(i, 2);
                [doa_candidate, localMin] = obj.golden_section(a, b, 20);
                if localMin
                    doa_candidates = [doa_candidates; doa_candidate];
                    mags = [mags; abs(obj.D(doa_candidate))];
                    idx = idx + 1;
                end
            end
            doa_angles = zeros(K, 1);
            for i = 1:K
                [~, idx] = min(mags);
                doa_angles(i) = doa_candidates(idx);
                mags = [mags(1:idx-1); mags(idx+1:end)];
                doa_candidates = [doa_candidates(1:idx-1); doa_candidates(idx+1:end)];
            end
        end

        function spec = D2(obj, angles)
            spec = zeros(1, length(angles));
            for i = 1:length(angles)
                spec(i) = obj.D(angles(i));
            end
        end

        function out = D(obj, theta)
            M = length(obj.gamma);
            z = exp(1i * pi * cosd(theta));
            out = 0;
            for i = 1:M
                for j = 1:M
                    out = out + obj.G(i, j) * z^(obj.gamma(j)-obj.gamma(i));
                end
            end
        end

        % Golden Section Search
        function [out, localMin] = golden_section(obj, a, b, ITER)
            g = 0.382;
            while true
                l = b - a;

                if l < 0.0005
                    break
                end

                a1 = a + g*l;
                b1 = b - g*l;
                if abs(obj.D(a1)) < abs(obj.D(b1))
                    b = b1;
                else
                    a = a1;
                end
            end
            out = 0.5 * (a + b);
            cond1 = (obj.D(out - 0.001) > obj.D(out));
            cond2 = (obj.D(out + 0.001) > obj.D(out));
            localMin = cond1 & cond2;
        end
    end
end