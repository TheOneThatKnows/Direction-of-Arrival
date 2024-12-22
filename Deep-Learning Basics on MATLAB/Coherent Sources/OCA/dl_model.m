%% CNN Structure

lgraph = [
    imageInputLayer([N N 2],"Normalization","none")
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
    fullyConnectedLayer(2*N-1)];
