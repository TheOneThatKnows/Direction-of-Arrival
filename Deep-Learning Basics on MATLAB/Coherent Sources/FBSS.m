function R_fb = FBSS(R, L)
M = size(R, 1);
P = M - L + 1;

J = zeros(M);
for i = 1:M
    J(i, M-i+1) = 1;
end
J

R_tilda = 0.5 * (R + J * conj(R) * J);

R_fb = zeros(L);
Z_eye = eye(M);
for p = 1:P
    Z = Z_eye(:, p:p+L-1)
    R_fb = R_fb + Z.' * R_tilda * Z;
end
R_fb = (1 / P) * R_fb;
end
