%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

M = 5; N = 10; Q = 121;

lgraph = [
    imageInputLayer([N N 3], "Normalization", "none")
    convolution2dLayer([4 4], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([3 3], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 128, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    flattenLayer
    fullyConnectedLayer(1024)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(256)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(128)
    leakyReluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(Q)
    fullyConnectedLayer(Q)
    reluLayer
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder II')
    convolution2dLayer([48, 1], 24, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([40, 1], 20, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([32, 1], 16, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([16, 1], 8, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([8, 1], 4, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([3, 1], 1, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    SSSpatialToChannelBroadcastLayer('Spatial Dimension Remover II')
    regressionLayer
    ];