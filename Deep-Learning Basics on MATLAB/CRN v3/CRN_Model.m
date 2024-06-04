%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

tempLayers = [
    imageInputLayer([5 5 1], "Name", "imageinput", "Normalization", "none")
    convolution2dLayer([3 3],10,"Name","conv")
    reluLayer("Name","relu")
    convolution2dLayer([2 2],5,"Name","conv_1")
    reluLayer("Name","relu_1")
    flattenLayer("Name","flatten")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    imageInputLayer([5 5 1], "Name", "imageinput_1", "Normalization", "none")
    convolution2dLayer([3 3],10,"Name","conv_2")
    reluLayer("Name","relu_2")
    convolution2dLayer([2 2],5,"Name","conv_3")
    reluLayer("Name","relu_3")
    flattenLayer("Name","flatten_1")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    depthConcatenationLayer(2,"Name","depthcat")
    fullyConnectedLayer(128, "Name", "fully_connected")
    reluLayer("Name","relu_4")
    fullyConnectedLayer(64, "Name", "fully_connected_1")
    reluLayer("Name","relu_5")
    fullyConnectedLayer(32, "Name", "fully_connected_2")
    reluLayer("Name","relu_6")
    fullyConnectedLayer(5, "Name", "fully_connected_3")
    sigmoidLayer("Name","sigmoid")
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Connect Layer Branches

lgraph = connectLayers(lgraph,"flatten","depthcat/in1");
lgraph = connectLayers(lgraph,"flatten_1","depthcat/in2");