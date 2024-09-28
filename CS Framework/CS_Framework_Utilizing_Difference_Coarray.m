% article link:
% https://ieeexplore.ieee.org/abstract/document/7744507

classdef CS_Framework_Utilizing_Difference_Coarray
    properties
        z;
        A;
        sg;
        mu;
    end
    methods
        function obj = CS_Framework_Utilizing_Difference_Coarray(z, A, mu)
            obj.z = z;
            obj.A = A;
            [~, Kg] = size(A);
            obj.sg = 0.001 * ones(Kg, 1);

            if nargin == 2
                obj.mu = [1; 0.1];
                return
            end

            obj.mu = mu;
        end

        function cost = Cost_Function(obj, sg)
            cost = obj.mu(1) * norm(sg, 1) + obj.mu(2) * norm(obj.z - obj.A * sg);
        end

        function gradient = Cost_Gradient(obj, sg)
            first_term = sign(sg);
            second_term = (1/norm(obj.z - obj.A * sg)) * ((sg' * (obj.A' * obj.A))' - (obj.z' * obj.A)');

            gradient = obj.mu(1) * first_term + obj.mu(2) * second_term;
        end
        
        function obj = Steepest_Descent_DOA(obj)
            cost_vector = zeros(500, 1);
            cost_idx = 1;
            while true
                cost = obj.Cost_Function(obj.sg);
                d = -obj.Cost_Gradient(obj.sg);
                alpha = obj.Dichotomous_Search(d, 0, 1);
                obj.sg = obj.sg + alpha * d;
                
                cost_vector(cost_idx) = cost;
                if cost_vector(20) ~= 0
                    local_mean = mean(cost_vector(cost_idx-19:cost_idx));
                    if abs(cost - local_mean) < 0.01 || cost_idx == 500
                        % cost
                        break
                    end
                end
                cost_idx = cost_idx + 1;
            end
        end

        function alpha = Dichotomous_Search(obj, d, a, b, LAMBDA)
            if nargin == 4
                LAMBDA = 0.01;
            end

            x1 = obj.sg + (0.5 * (a + b) - LAMBDA) * d;
            x2 = obj.sg + (0.5 * (a + b) + LAMBDA) * d;

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

                x1 = obj.sg + (0.5 * (a + b) - LAMBDA) * d;
                x2 = obj.sg + (0.5 * (a + b) + LAMBDA) * d;
            end

            alpha = 0.5 * (a + b);
        end
    end
end