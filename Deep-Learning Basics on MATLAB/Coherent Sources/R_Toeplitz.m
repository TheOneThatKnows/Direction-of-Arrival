function R_out = R_Toeplitz(R, half_or_full)
M = size(R, 1);

if half_or_full == "half"
    count = (M-1 - rem(M-1, 2))/2 + 1;
else
    count = M;
end
R_ = zeros(M, M, count);
R_m = zeros(M);

C_ = zeros(M, M, count);
C_m = zeros(M);

for m = 0:count-1
    J = zeros(M - m);
    for j = 1:M-m
        J(j, M-m-j+1) = 1;
    end

    % R
    r_m = [zeros(1, m) (R(m+1:M, m+1).' * J) R(m+1, m+2:M) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        R_m(i+1, :) = r_m(start_idx:end_idx).';
    end

    R_(:, :, m+1) = R_m;

    % C
    c_m = [zeros(1, m) ones(1, 2*M-2*m-1) zeros(1, m)].';
    
    for i = 0:M-1
        start_idx = M - i;
        end_idx = 2 * M - i - 1;
        C_m(i+1, :) = c_m(start_idx:end_idx).';
    end

    C_(:, :, m+1) = C_m;
end

R_sum = zeros(M);
C_sum = zeros(M);
for m = 0:(M-1)/2
    R_sum = R_sum + R_(:, :, m+1);
    C_sum = C_sum + C_(:, :, m+1);
end

R_out = R_sum ./ C_sum;
end