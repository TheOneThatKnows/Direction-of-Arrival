classdef CS_Framework
    properties
        y;
        A;
        Sg;
        mu;
        d;
    end
    methods
        function obj = CS_Framework(y, A)
            obj.y = y;
            obj.A = A;
            [~, P] = size(y);
            [~, Kg] = size(A);
            obj.Sg = zeros(Kg, P);
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
                sg_hat(i) = norm(Sg(i, :));
            end
        end

        function gradient = Cost_Gradient(obj, Sg, sg_hat)
            [Kg, P] = size(Sg);
            
            first_term = zeros(Kg, P);
            for i = 1:P
                for j = 1:Kg
                    first_term(j, i) = Sg(j, i) / sg_hat(j);
                end
            end

            second_term = zeros(Kg, P);
            V = obj.y - obj.A * Sg;
            norm_value = norm(V(:));
            for i = 1:P
                for j = 1:Kg
                    second_term(j, i) = -(1/norm_value) * (obj.A(:, j).' * V(:, i));
                end
            end

            % first_term = zeros(P*Kg, 1);
            % for i = 1:P
            %     for j = 1:Kg
            %         idx = (i-1) * Kg + j;
            %         first_term(idx) = Sg(j, i) / sg_hat(j);
            %     end
            % end
            % 
            % second_term = zeros(P*Kg, 1);
            % V = obj.y - obj.A * Sg;
            % norm_value = norm(V(:));
            % for i = 1:P
            %     for j = 1:Kg
            %         idx = (i-1) * Kg + j;
            %         second_term(idx) = -(1/norm_value) * (obj.A(:, j)' * V(:, i));
            %     end
            % end

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
                obj.d = -obj.Cost_Gradient(obj.Sg, sg_hat);
                alpha = obj.Dichotomous_Search(0, 1);
                obj.Sg = obj.Sg + alpha * obj.d;
                
                if abs(cost - old_cost) < 0.01
                    break
                end

                old_cost = cost;
            end
        end

        function alpha = Dichotomous_Search(obj, a, b, LAMBDA)
            if nargin == 3
                LAMBDA = 0.01;
            end

            x1 = obj.Sg + (0.5 * (a + b) - LAMBDA) * obj.d;
            x2 = obj.Sg + (0.5 * (a + b) + LAMBDA) * obj.d;

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

                x1 = obj.Sg + (0.5 * (a + b) - LAMBDA) * obj.d;
                x2 = obj.Sg + (0.5 * (a + b) + LAMBDA) * obj.d;
            end

            alpha = 0.5 * (a + b);
        end
    end
end