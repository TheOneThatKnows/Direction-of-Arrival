%% Create Layer Graph

M = 5; Q = 121;

lgraph = [
    featureInputLayer((M-1) * (M-1), "Normalization", "none")
    fullyConnectedLayer((M-2) * (M-2))
    fullyConnectedLayer((M-1) * (M-1))
    regressionLayer
    ];