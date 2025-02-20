classdef LRMD_Object
    properties
        A;
        B;
        W;
        R;
        gamma;
        beta;
    end
    methods
        function obj = LRMD_Object(A, B, W, R, gamma, beta)
            obj.A = A;
            obj.B = B;
            obj.W = W;
            obj.R = R;
            obj.gamma = gamma;
            obj.beta = beta;
        end

        function cost = Cost(obj, A, B)
            h = 0.5 * norm(obj.W .* (obj.R - A * B), "fro")^2;
            f = 0;
            for j = 1:size(A, 2)
                f = f + norm(A(:, j));
            end
            f = obj.gamma * f;
            g = 0.5 * obj.gamma * obj.beta * norm(B, "fro")^2;
            cost = h + f + g;
        end

        function [A_k1, B_k1] = Steepest_Descent(obj, ITER)
            N = size(obj.A, 1);
            d = size(obj.A, 2);

            grad_A = zeros(N, d);
            grad_B = zeros(d, N);

            A_k1 = obj.A;
            B_k1 = obj.B;

            for k = 1:ITER
                A_k = A_k1;
                B_k = B_k1;

                % Update A
                for i = 1:N
                    for j = 1:d
                        grad_A(i, j) = -0.5 * B_k(j, :) * (obj.W(i, :) .* (obj.R(i, :) - A_k(i, :) * B_k))';
                    end
                end
                alpha = obj.Dichotomous_Search_A(-grad_A, 0, 100);
                A_k1 = A_k - alpha * grad_A;

                % Update B
                for i = 1:d
                    for j = 1:N
                        grad_B(i, j) = -0.5 * (obj.W(:, j) .* (obj.R(:, j) - A_k1 * B_k(:, j)))' * A_k1(:, i);
                    end
                end
                alpha = obj.Dichotomous_Search_B(-grad_B, 0, 100);
                B_k1 = B_k - alpha * grad_B;
            end
        end

        function alpha = Dichotomous_Search_A(obj, d, a, b, LAMBDA)
            if nargin == 4
                LAMBDA = 0.01;
            end

            x1 = obj.A + (0.5 * (a + b) - LAMBDA) * d;
            x2 = obj.A + (0.5 * (a + b) + LAMBDA) * d;

            while true
                cost1 = obj.Cost(x1, obj.B);
                cost2 = obj.Cost(x2, obj.B);

                if cost1 > cost2
                    a = (0.5 * (a + b) - LAMBDA);
                else
                    b = (0.5 * (a + b) + LAMBDA);
                end

                if b - a < 0.05
                    break
                end

                x1 = obj.A + (0.5 * (a + b) - LAMBDA) * d;
                x2 = obj.A + (0.5 * (a + b) + LAMBDA) * d;
            end

            alpha = 0.5 * (a + b);
        end
        function alpha = Dichotomous_Search_B(obj, d, a, b, LAMBDA)
            if nargin == 4
                LAMBDA = 0.01;
            end

            x1 = obj.B + (0.5 * (a + b) - LAMBDA) * d;
            x2 = obj.B + (0.5 * (a + b) + LAMBDA) * d;

            while true
                cost1 = obj.Cost(obj.A, x1);
                cost2 = obj.Cost(obj.A, x2);

                if cost1 > cost2
                    a = (0.5 * (a + b) - LAMBDA);
                else
                    b = (0.5 * (a + b) + LAMBDA);
                end

                if b - a < 0.05
                    break
                end

                x1 = obj.B + (0.5 * (a + b) - LAMBDA) * d;
                x2 = obj.B + (0.5 * (a + b) + LAMBDA) * d;
            end

            alpha = 0.5 * (a + b);
        end
    end
end