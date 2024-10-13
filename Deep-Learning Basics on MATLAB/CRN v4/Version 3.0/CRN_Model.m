%% Add Layer Branches

M = 5; N = 10;

lgraph = [
    imageInputLayer([M M 3], "Normalization", "none")
    convolution2dLayer([2 2], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    flattenLayer
    fullyConnectedLayer(512)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(256)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(128)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(2*N-1)
    % leakyReluLayer
    regressionLayer
    ];