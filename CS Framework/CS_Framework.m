% article link:
% https://ieeexplore.ieee.org/abstract/document/7744507

classdef CS_Framework
    properties
        y;
        A;
        Sg;
        mu;
    end
    methods
        function obj = CS_Framework(y, A)
            obj.y = y;
            obj.A = A;
            [~, P] = size(y);
            [~, Kg] = size(A);
            obj.Sg = 0.001 * ones(Kg, P);
        end

        function [cost, sg_hat] = Cost_Function(obj, Sg)
            if nargin == 1
                Sg = obj.Sg;
            end
            sg_hat = obj.Calculate_sg_hat(Sg);
            vectorized_term = obj.y - obj.A * Sg;
            vectorized_term = vectorized_term(:);
            cost = obj.mu(1) * norm(sg_hat, 1) + obj.mu(2) * norm(vectorized_term);
        end

        function sg_hat = Calculate_sg_hat(~, Sg)
            [Kg, ~] = size(Sg);
            sg_hat = zeros(Kg, 1);
            for i = 1:Kg
                sg_hat(i) = norm(Sg(i, :))^2;
            end
        end

        function gradient = Cost_Gradient(obj, Sg, sg_hat)
            [Kg, P] = size(Sg);
            
            first_term = zeros(Kg, P);
            for i = 1:P
                for j = 1:Kg
                    % first_term(j, i) = 0.5 * conj(Sg(j, i)) / sg_hat(j);
                    first_term(j, i) = conj(Sg(j, i));
                end
            end

            second_term = zeros(Kg, P);
            V = obj.y - obj.A * Sg;
            norm_value = norm(V(:));
            for i = 1:P
                for j = 1:Kg
                    second_term(j, i) = -0.5 * (1/norm_value) * (V(:, i)' * obj.A(:, j));
                end
            end

            gradient = obj.mu(1) * first_term + obj.mu(2) * second_term;
        end
        
        function obj = Steepest_Descent_DOA(obj, mu)
            if nargin == 1
                obj.mu = [1; 0.1];
            else
                obj.mu = mu;
            end

            old_cost = inf;

            while true
                [cost, sg_hat] = obj.Cost_Function();
                d = -obj.Cost_Gradient(obj.Sg, sg_hat);
                alpha = obj.Dichotomous_Search(d, 0, 1);
                obj.Sg = obj.Sg + alpha * d;
                
                if abs(cost - old_cost) < 0.005
                    break
                end

                old_cost = cost;
            end
        end

        function alpha = Dichotomous_Search(obj, d, a, b, LAMBDA)
            if nargin == 4
                LAMBDA = 0.01;
            end

            x1 = obj.Sg + (0.5 * (a + b) - LAMBDA) * d;
            x2 = obj.Sg + (0.5 * (a + b) + LAMBDA) * d;

            while true
                f1 = obj.Cost_Function(x1);
                f2 = obj.Cost_Function(x2);

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
        end
    end
end