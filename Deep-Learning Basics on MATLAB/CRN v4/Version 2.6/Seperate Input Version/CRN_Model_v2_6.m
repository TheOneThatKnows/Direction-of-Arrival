%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

M = 5; N = 10; Q = 121;

tempLayers = [
    imageInputLayer([M M 1], "Normalization", "none")
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    flattenLayer("Name","flatten")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    imageInputLayer([M M 1], "Normalization", "none")
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    flattenLayer("Name","flatten_1")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    imageInputLayer([M M 1], "Normalization", "none")
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    convolution2dLayer([2 2], 32, "Stride", 1, "Padding", 0)
    batchNormalizationLayer
    leakyReluLayer
    flattenLayer("Name","flatten_2")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    depthConcatenationLayer(3, "Name","depthcat")
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
    regressionLayer];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Connect Layer Branches

lgraph = connectLayers(lgraph,"flatten","depthcat/in1");
lgraph = connectLayers(lgraph,"flatten_1","depthcat/in2");
lgraph = connectLayers(lgraph,"flatten_2","depthcat/in3");

%% Plot Layers

figure;
plot(lgraph);