%% Create Layer Graph

M = 5; N = 10;

lgraph = [
    featureInputLayer(2 * M, "Normalization", "none")
    fullyConnectedLayer(512)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(256)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(128)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(64)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(32)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(2 * N)
    regressionLayer
    ];