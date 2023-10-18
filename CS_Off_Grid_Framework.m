% article link:
% https://www.sciencedirect.com/science/article/pii/B9780128118870000110

classdef CS_Off_Grid_Framework
    properties
        y;
        A;
        B;
        Sg;
        beta;
        mu;
    end
    methods
        function obj = CS_Off_Grid_Framework(y, A, B)
            obj.y = y;
            obj.A = A;
            obj.B = B;
            [~, P] = size(y);
            [~, Kg] = size(A);
            obj.Sg = 0.001 * ones(Kg, P);
            obj.beta = 0.0001 * ones(Kg, 1);
        end

        function cost = Cost_Function(obj, Sg, beta)
            sg_hat = obj.Calculate_sg_hat(Sg);
            Phi = obj.A + obj.B * diag(beta);
            vectorized_term = obj.y - Phi * Sg;
            vectorized_term = vectorized_term(:);
            cost = obj.mu(1) * sum(sg_hat) + 0.5 * (vectorized_term' * vectorized_term) + obj.mu(2) * (beta' * beta);
        end

        function sg_hat = Calculate_sg_hat(~, Sg)
            [Kg, ~] = size(Sg);
            sg_hat = zeros(Kg, 1);
            for i = 1:Kg
                sg_hat(i) = norm(Sg(i, :))^2;
            end
        end

        function [gradient_Sg, gradient_beta] = Cost_Gradient(obj, Sg)
            [Kg, P] = size(Sg);
            
            first_term = zeros(Kg, P);
            for i = 1:P
                for j = 1:Kg
                    first_term(j, i) = conj(Sg(j, i));
                end
            end

            second_term = zeros(Kg, P);
            Phi = obj.A + obj.B * diag(obj.beta);
            V = obj.y - Phi * Sg;
            % norm_value = norm(V(:));
            for i = 1:P
                for j = 1:Kg
                    second_term(j, i) = -0.5 * (V(:, i)' * Phi(:, j));
                end
            end

            gradient_Sg = obj.mu(1) * first_term + second_term;

            first_term = zeros(Kg, 1);
            for i = 1:Kg
                first_term(i) = -0.5 * Sg(i, :) * V' * obj.B(:, i);
            end
            second_term = conj(obj.beta);

            gradient_beta = first_term + obj.mu(2) * second_term;
        end
        
        function obj = Steepest_Descent_DOA(obj, mu)
            if nargin == 1
                obj.mu = [1; 0.1];
            else
                obj.mu = mu;
            end

            old_cost = inf;

            while true
                while true
                    cost = obj.Cost_Function(obj.Sg, obj.beta);
                    [g, ~] = obj.Cost_Gradient(obj.Sg);
                    d = -g;
                    alpha = obj.Dichotomous_Search("Sg", d, 0, 1);
                    obj.Sg = obj.Sg + alpha * d;
                    cost
                    if abs(cost - old_cost) < 0.005
                        break
                    end

                    old_cost = cost;
                end
                [~, g] = obj.Cost_Gradient(obj.Sg);
                d = -g;
                alpha = obj.Dichotomous_Search("beta", d, 0, 1);
                obj.beta = obj.beta + alpha * d;

                if abs(cost - old_cost) < 0.005
                    break
                end

                old_cost = cost;
            end
        end

        function alpha = Dichotomous_Search(obj, mode, d, a, b, LAMBDA)
            if nargin == 5
                LAMBDA = 0.01;
            end

            if mode == "Sg"
                x1 = obj.Sg + (0.5 * (a + b) - LAMBDA) * d;
                x2 = obj.Sg + (0.5 * (a + b) + LAMBDA) * d;

                while true
                    f1 = obj.Cost_Function(x1, obj.beta);
                    f2 = obj.Cost_Function(x2, obj.beta);

                    if f1 > f2
                        a = (0.5 * (a + b) - LAMBDA);
                    else
                        b = (0.5 * (a + b) + LAMBDA);
                    end

                    if b - a < 0.05
                        break
                    end

                    x1 = obj.Sg + (0.5 * (a + b) - LAMBDA) * d;
                    x2 = obj.Sg + (0.5 * (a + b) + LAMBDA) * d;
                end

                alpha = 0.5 * (a + b);
            else % if mode == "beta"
                x1 = obj.beta + (0.5 * (a + b) - LAMBDA) * d;
                x2 = obj.beta + (0.5 * (a + b) + LAMBDA) * d;

                while true
                    f1 = obj.Cost_Function(obj.Sg, x1);
                    f2 = obj.Cost_Function(obj.Sg, x2);

                    if f1 > f2
                        a = (0.5 * (a + b) - LAMBDA);
                    else
                        b = (0.5 * (a + b) + LAMBDA);
                    end

                    if b - a < 0.05
                        break
                    end

                    x1 = obj.beta + (0.5 * (a + b) - LAMBDA) * d;
                    x2 = obj.beta + (0.5 * (a + b) + LAMBDA) * d;
                end

                alpha = 0.5 * (a + b);
            end
        end
    end
end