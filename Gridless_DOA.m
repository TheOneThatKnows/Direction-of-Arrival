% article link:
% https://ieeexplore.ieee.org/abstract/document/9384289

classdef Gridless_DOA
    properties
        G;
        gamma;
    end
    methods
        function obj = Gridless_DOA()
        end

        % Irregular Vandermonde Decomposition
        function [obj, doa_angles, z, c] = IVD(obj, T, gamma, K)
            % T: Autocovarience Matrix
            % gamma: Sensor Positions
            % K: # of sources

            obj.gamma = gamma;
            M = length(gamma);  % # of sensors
            [U, ~, ~] = svd(T);
            U_N = U(:, K+1:M);  % Noise Space

            [obj, doa_angles] = obj.estimate_doas(U_N, M, K);
            z = exp(1i * pi * cosd(doa_angles));
            W = obj.IVM(gamma, z);
            W_pseudo_inv = (W' * W)\W';
            C = W_pseudo_inv * T * W_pseudo_inv';
            c = zeros(K, 1);
            for i = 1:K
                c(i) = C(i, i);
            end
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
                [doa_candidate, localMin] = obj.golden_section(a, b);
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
        function [out, localMin] = golden_section(obj, a, b)
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

        % Irregular Vandermonde Matrix
        function W = IVM(~, gamma, z)
            M = length(gamma);
            K = length(z);

            W = zeros(M, K);
            for i = 1:K
                W(:, i) = z(i).^(gamma.');
            end
        end

        % Projection Onto Irregular Toeplitz Set
        function [obj, T_projected] = PTG(obj, T, gamma, K)
            [obj, ~, z, c] = obj.IVD(T, gamma, K);
            W = obj.IVM(gamma, z);
            T_projected = W * diag(c) * W';
        end

        % Projection Onto Positive Semi-Definite Cone
        function M_out = PPSD(~, M_in)
            [N, ~] = size(M_in);
            M_out = zeros(N);
            [eig_vec, eig_val] = eig(M_in, 'vector');
            
            for i = 1:N
                M_out = M_out + max(0, eig_val(i)) * (eig_vec(:, i) * eig_vec(:, i)');
            end
        end

        % Projection Onto S(T_gamma_K, Y)
        function M_out = PS(obj, M_in, Y, gamma, K)
            M_out = zeros(size(M_in));
            [rY, ~] = size(Y);
            B = M_in(1:rY, 1:rY);
            Q = M_in(rY+1:end, rY+1:end);

            M_out(1:rY, 1:rY) = obj.PTG(B, gamma, K);
            M_out(1:rY, rY+1:end) = Y;
            M_out(rY+1:end, 1:rY) = Y';
            M_out(rY+1:end, rY+1:end) = Q;
        end

        % AP Gridless
        function [obj, doa_angles] = AP_Gridless(obj, Y, gamma, K)
            [M, N] = size(Y);
            B = zeros(M);
            Q = (1 / N) * (Y' * Y);
            L = [B Y; Y' Q];
            while true
                H = obj.PPSD(L);
                T = H(1:M, 1:M);
                L_old = L;
                [obj, T_projected] = obj.PTG(T, gamma, K);
                L = [T_projected Y; Y' Q];

                if norm(L - L_old) <= 1e-7
                    break
                end
            end
            T = L(1:M, 1:M);
            [obj, doa_angles, ~, ~] = obj.IVD(T, gamma, K);
        end
    end
end