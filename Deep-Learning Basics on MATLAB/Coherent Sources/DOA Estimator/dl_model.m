%% CNN Structure

lgraph = [
    imageInputLayer([N N 2],"Normalization","none")
    convolution2dLayer([2 2],128)
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([2 2],128)
    batchNormalizationLayer
    reluLayer
    convolution2dLayer([2 2],128)
    batchNormalizationLayer
    reluLayer
    flattenLayer
    fullyConnectedLayer(1024)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(512)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(128)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(Q)];