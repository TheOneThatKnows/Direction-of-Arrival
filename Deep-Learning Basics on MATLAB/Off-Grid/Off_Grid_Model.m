%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

M = 5; Q = 121;

lgraph = [
    featureInputLayer((M-2) * (M-2), "Normalization", "none")
    fullyConnectedLayer(16 * (M-2))
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(32 * (M-2))
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(64 * (M-2))
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(128 * (M-2))
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.3)
    fullyConnectedLayer(128 * (M-2))
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(Q)
    reluLayer
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder')
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
    SSSpatialToChannelBroadcastLayer('Spatial Dimension Remover')
    regressionLayer
    ];