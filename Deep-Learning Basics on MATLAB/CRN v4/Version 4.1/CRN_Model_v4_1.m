%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

lgraph = layerGraph();

%% Create Layer Graph

M = 5; N = 10; Q = 121;

tempLayers = [
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
    fullyConnectedLayer(Q)
    batchNormalizationLayer
    reluLayer
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder I')
    convolution2dLayer([21, 1], 24, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([15, 1], 20, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([11, 1], 12, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([5, 1], 5, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([3, 1], 1, "Padding", "same")
    batchNormalizationLayer
    reluLayer("Name", "concat_in_1")
    ];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    featureInputLayer(Q, "Normalization", "none")
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder II')
    ];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    concatenationLayer(3,2,"Name","concat")
    convolution2dLayer([21, 1], 24, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([15, 1], 20, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([11, 1], 12, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([5, 1], 5, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([3, 1], 1, "Padding", "same")
    batchNormalizationLayer
    reluLayer
    SSSpatialToChannelBroadcastLayer('Spatial Dimension Remover')
    regressionLayer
    ];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Connect Layer Branches

lgraph = connectLayers(lgraph,"concat_in_1","concat/in1");
lgraph = connectLayers(lgraph,"Spatial Dimension Adder II","concat/in2");

%% Plot Layers

figure;
plot(lgraph);