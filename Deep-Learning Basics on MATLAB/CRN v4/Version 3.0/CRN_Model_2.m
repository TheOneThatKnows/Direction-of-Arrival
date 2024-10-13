%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

N = 10; Q = 121;

layers = [
    featureInputLayer(2 * N - 1, "Normalization", "none")
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

lgraph = addLayers(lgraph, layers);

clear layers

%% Plot Layers

figure;
plot(lgraph);