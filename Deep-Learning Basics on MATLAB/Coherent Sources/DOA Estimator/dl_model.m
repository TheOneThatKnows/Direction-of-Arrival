%% Custom Layers

addpath(['D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival' ...
    '\Deep-Learning Basics on MATLAB\Custom Layers']);

%% CNN Structure

lgraph = [
    imageInputLayer([N N 3],"Normalization","none")
    convolution2dLayer([4 4],128)
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([3 3],128)
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([3 3],128)
    batchNormalizationLayer
    reluLayer
    flattenLayer
    fullyConnectedLayer(656)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(374)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(213)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(Q)
    reluLayer
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder II')
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
    SSSpatialToChannelBroadcastLayer('Spatial Dimension Remover II')];
