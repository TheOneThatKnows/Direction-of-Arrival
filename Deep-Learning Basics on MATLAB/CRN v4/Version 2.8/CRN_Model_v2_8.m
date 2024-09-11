%% Add Custom Layer Library
% conv1D can be discarded, combinations can be applied (1D->MLP, MLP->1D, ...)
addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

M = 5; N = 10; Q = 121;

tempLayers = [
    imageInputLayer([M M 3], "Name", "imageinput", "Normalization", "none")
    convolution2dLayer([3 3],10,"Name","conv")
    reluLayer("Name","relu")
    convolution2dLayer([2 2],5,"Name","conv_1")
    reluLayer("Name","relu_1")
    flattenLayer("Name","flatten")
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(Q)
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Not used

tempLayers = [
    depthConcatenationLayer(2,"Name","depthcat")
    fullyConnectedLayer(40, "Name", "fully_connected")
    reluLayer("Name","relu_4")
    SpatialToChannelBroadcastLayer('Spatial Dimension Adder')
    convolution2dLayer([15 1],3,"Name","conv1d","Padding","same")
    reluLayer("Name","relu_5")
    convolution2dLayer([5 1],1,"Name","conv1d_1","Padding","same")
    reluLayer("Name","relu_6")
    SSSpatialToChannelBroadcastLayer('Spatial Dimension Remover')
    fullyConnectedLayer(2 * N - 1, "Name","fully_connected_2")
    % reluLayer("Name","relu_7")
    % sigmoidLayer("Name","sigmoid")
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(Q)
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;