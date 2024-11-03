%% LSTM Structure

numChannels = size(data{1}, 2);
numHiddenUnits = [125 120];
numClasses = N - 1;

layers = [
    sequenceInputLayer(numChannels)
    bilstmLayer(numHiddenUnits(1),OutputMode="sequence")
    dropoutLayer(0.2)
    bilstmLayer(numHiddenUnits(2),OutputMode="last")
    dropoutLayer(0.2)
    fullyConnectedLayer(numClasses)
    softmaxLayer];

%% CNN Structure

lgraph = [
    imageInputLayer([M M 2],"Normalization","none")
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
    fullyConnectedLayer(N-1)
    softmaxLayer];