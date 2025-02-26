from numpy import concatenate, argsort, sort, arange
from numpy import zeros, ones, eye
from numpy import pi, inf
from numpy import abs, sign, mean, var, sqrt, cos, exp, conj, dot, round
from numpy.random import randn, rand, random, permutation
from numpy.linalg import eig, inv
from scipy.linalg import toeplitz
from scipy.signal import find_peaks

# DOA Estimator
def DOA_Estimate(spec, angle_spec, K):
    # Extend the spectrum with -inf at both ends
    spec = concatenate(([-inf], spec, [-inf]))

    # Find peaks in the spectrum
    peak_inds, _ = find_peaks(spec)
    peak_mags = spec[peak_inds]

    # Initialize the DOA estimate array
    doa_est = zeros(K)

    # Sort the peaks by magnitude in descending order
    sorted_inds = argsort(peak_mags)[::-1]
    peak_inds = peak_inds[sorted_inds]

    # Ensure peak_inds has enough elements
    peak_inds = concatenate((peak_inds, zeros(K - 1)))

    # First DOA estimate
    idx = peak_inds.astype(int)
    doa_est = angle_spec[idx - 1]

    # Loop to estimate remaining DOAs
    for i in range(1, K):
        idx = int(peak_inds[i])
        if idx == 0:
            doa_est[i:K] = mean(doa_est[:i])
            break
        else:
            doa_est[i] = angle_spec[idx - 1]

    # Sort the DOA estimates
    doa_est = doa_est[:K]
    doa_est = sort(doa_est)

    return doa_est

def DOA_Estimate_Local(spec, angle_spec, doa, local_delta):
    delta_spec = angle_spec[1] - angle_spec[0]
    doa_inds = round((doa - angle_spec[0]) / delta_spec)
    local_delta = round(local_delta / delta_spec) + 1

    # Initialize the DOA estimate array
    K = len(doa)
    doa_est = zeros(K)

    for i in range(K):
        ind_min = max(0, int(doa_inds[i]-local_delta))
        ind_max = min(len(spec), int(doa_inds[i]+local_delta))

        doa_est[i] = DOA_Estimate(spec[ind_min:ind_max], angle_spec[ind_min:ind_max], 1)[0]

    # Sort the DOA estimates
    doa_est = sort(doa_est)

    return doa_est

def DOA_Generate(K, phi_min, phi_max, delta_phi):
    doa = (-2 * delta_phi) * ones(K)
    i = 0
    while True:
        temp_angle = phi_min + random() * (phi_max - phi_min)
        temp_array = abs(doa - temp_angle)
        if any(temp_array < delta_phi * 2):
            i -= 1
        else:
            doa[i] = temp_angle
        if i is K-1:
            break
        i += 1

    return sort(doa)

def Array_Manifold(sensor_locations, doa):
    M = sensor_locations.shape[1]
    K = len(doa)

    A = zeros((M, K), dtype=complex)
    for i in range(K):
        A[:, i] = exp(1j * pi * sensor_locations.T * cos(doa[i])).squeeze()

    return A

def Source_Generate(K, K_coherent, L, vars, alpha=None):
    if alpha is None:
        alpha = concatenate(([[1]], (0.75 + 0.25 * rand((sign(K_coherent) * (K_coherent - 1)), 1)) *
                                 exp(1j * pi * rand((sign(K_coherent) * (K_coherent - 1)), 1))))

    v = zeros((K, L), dtype=complex)
    v[:K_coherent, :] = ones((K_coherent, 1)) @ (randn(1, L) + 1j * randn(1, L))
    v[K_coherent:, :] = randn(K - K_coherent, L) + 1j * randn(K - K_coherent, L)

    s = v * exp(1j * (2 * pi * 0.5 * arange(L)))

    for i in range(K):
        s[i, :] /= sqrt(var(s[i, :]))

    s = concatenate((sqrt(vars[0]) * alpha, sqrt(vars[1:])), axis=0) * s

    shuffled_inds = permutation(K)
    s = s[shuffled_inds, :]

    return s

def Noise_Generate(M, L, SNR_dB):
    return (1 / sqrt(2 * pow(10, SNR_dB / 10))) * (randn(M, L) + 1j * randn(M, L))

def MUSIC(K, R, sensor_locations, angles):
    spec = zeros(len(angles))

    M_v = R.shape[0]

    # Eigen decomposition of the covariance matrix
    eig_vals, eig_vecs = eig(R)
    sorted_inds = argsort(eig_vals)
    eig_vecs = eig_vecs[:, sorted_inds]
    G = eig_vecs[:, :M_v-K]  # noise space

    for i in range(len(angles)):
        a_ = exp(1j * pi * sensor_locations.T * cos(angles[i]))
        spec[i] = (1 / abs(conj(a_).T @ (G @ conj(G).T) @ a_)).item()

    spec = spec / max(spec)
    return spec

def CBF(R, sensor_locations, angles):
    spec = zeros(len(angles))
    M = sensor_locations.shape[1]

    for i in range(len(angles)):
        h = exp(1j * pi * sensor_locations.T * cos(angles[i]))
        spec[i] = abs(conj(h).T @ R @ h).item()

    spec = spec / max(spec)
    return spec

def Capon(R, sensor_locations, angles):
    spec = zeros(len(angles))
    M = sensor_locations.shape[1]

    for i in range(len(angles)):
        a = exp(1j * pi * sensor_locations.T * cos(angles[i]))
        h = (inv(R) @ a) / (conj(a).T @ inv(R) @ a)
        spec[i] = abs(conj(h).T @ R @ h).item()

    spec = spec / max(spec)
    return spec

def Sensor_Placement(sensor_locations):
    sensor_placement = zeros((1, sensor_locations[0, -1] + 1))
    sensor_placement[0, sensor_locations[0, :]] = 1

    return sensor_placement

def Diff_Coarray(sensor_locations):
    sensor_placement = Sensor_Placement(sensor_locations)
    N = sensor_placement.shape[1]
    diff_coarray = zeros((1, 2*N-1))
    for i in range(N):
        diff_coarray[0, N-1+i] = sensor_placement[0, i:] @ sensor_placement[0, :N-i].T
    diff_coarray[0, :N-1] = diff_coarray[0, arange(2*N-2, N-1, -1)]

    return diff_coarray

def Uniform_Degrees_Of_Freedom(sensor_locations, one_sided=False):
    diff_coarray = Diff_Coarray(sensor_locations)
    uniform_degrees_of_freedom = 0
    N = sensor_locations[0, -1] + 1
    for i in range(N):
        if diff_coarray[0, N-1+i] == 0:
            break
        uniform_degrees_of_freedom += 1
    if not one_sided:
        uniform_degrees_of_freedom = 2 * uniform_degrees_of_freedom - 1

    return uniform_degrees_of_freedom

def R_Toeplitz(R, half_or_full):
    M = R.shape[0]

    if half_or_full == "half":
        count = (M - 1 - (M - 1) % 2) // 2 + 1
    else:
        count = M

    R_ = zeros((count, M, M), dtype=complex)
    R_m = zeros((M, M), dtype=complex)

    C_ = zeros((count, M, M), dtype=complex)
    C_m = zeros((M, M), dtype=complex)

    for m in range(count):
        J = zeros((M - m, M - m))
        for j in range(M - m):
            J[j, M - m - j - 1] = 1

        # R
        r_m = concatenate((zeros(m), (R[m:M, m] @ J).flatten(), R[m, m + 1:M], zeros(m)))
        
        for i in range(M):
            start_idx = M - i - 1
            end_idx = 2 * M - i - 1
            R_m[i, :] = r_m[start_idx:end_idx]

        R_[m, :, :] = R_m

        # C
        c_m = concatenate((zeros(m), ones(2 * M - 2 * m - 1), zeros(m)))
        
        for i in range(M):
            start_idx = M - i - 1
            end_idx = 2 * M - i - 1
            C_m[i, :] = c_m[start_idx:end_idx]

        C_[m, :, :] = C_m

    R_sum = zeros((M, M), dtype=complex)
    C_sum = zeros((M, M), dtype=complex)
    for m in range(count):
        R_sum += R_[m, :, :]
        C_sum += C_[m, :, :]

    R_out = R_sum / C_sum
    return R_out

def Virtual_Covariance_Column(R, sensor_locations, outMatrix=False):
    M = sensor_locations.shape[1]
    N = sensor_locations[0, -1] + 1
    diff_coarray = Diff_Coarray(sensor_locations)

    sparse_entries = zeros((2*N-1, 1), dtype=complex)
    for i in range(M):
        for j in range(M):
            diff = sensor_locations[0, i] - sensor_locations[0, j]
            count = diff_coarray[0, N - 1 + diff]
            sparse_entries[N - 1 + diff, 0] = sparse_entries[N - 1 + diff, 0] + (1 / count) * R[i, j]

    column = concatenate((sparse_entries[N-1:N], 0.5 * (sparse_entries[N:] + conj(sparse_entries[arange(N-2, -1, -1)]))))

    if outMatrix:
        return toeplitz(column)
    return column

def Spatial_Smoothing(R, sensor_locations):
    M = sensor_locations.shape[1]
    N = sensor_locations[0, -1] + 1
    diff_coarray = Diff_Coarray(sensor_locations)
    uDOF = Uniform_Degrees_Of_Freedom(sensor_locations, True)

    sparse_entries = zeros((2*uDOF-1, 1), dtype=complex)
    for i in range(M):
        for j in range(M):
            diff = sensor_locations[0, i] - sensor_locations[0, j]
            if abs(diff) > uDOF - 1:
                continue
            count = diff_coarray[0, N - 1 + diff]
            sparse_entries[uDOF - 1 + diff, 0] = sparse_entries[uDOF - 1 + diff, 0] + (1 / count) * R[i, j]
    
    R_out = zeros((uDOF, uDOF),dtype=complex)
    for i in range(uDOF):
        z1_i = sparse_entries[i:uDOF + i]
        R_out = R_out + (1 / uDOF) * (z1_i @ conj(z1_i).T)

    return R_out

def Vector2Matrix(r, sensor_locations):
    M = sensor_locations.shape[1]

    R_out = zeros((M, M), dtype=complex)
    for i in range(M):
        R_out[i, i] = r[0]
    R_out /= 2
    for i in arange(1, M):
        for j in range(i):
            diff = sensor_locations[0, i] - sensor_locations[0, j]
            R_out[i, j] = r[diff]
    R_out += conj(R_out).T

    return R_out

def FBSS(R, subarray_size):
    M = R.shape[0]
    subarray_count = M - subarray_size + 1

    J = zeros((M, M))
    for i in range(M):
        J[i, M-i-1] = 1
    
    R_tilda = 0.5 * (R + J @ R.T @ J)

    R_fb = zeros((subarray_size, subarray_size), dtype=complex)
    Z_eye = eye(M)
    for p in range(subarray_count):
        Z = Z_eye[:, p:p+subarray_size]
        R_fb += Z.T @ R_tilda @ Z
        
    R_fb = (1 / subarray_count) * R_fb
    return R_fb