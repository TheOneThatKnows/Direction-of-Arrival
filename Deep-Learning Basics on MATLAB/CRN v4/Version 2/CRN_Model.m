%% Add Custom Layer Library

addpath('D:\D\Alp\Master ODTÜ\Thesis\DOA\Codes\Direction-of-Arrival\Deep-Learning Basics on MATLAB\Custom Layers');

%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

M = 5; N = 10;

tempLayers = [
    imageInputLayer([M M 1], "Name", "imageinput", "Normalization", "none")
    convolution2dLayer([3 3],10,"Name","conv")
    reluLayer("Name","relu")
    convolution2dLayer([2 2],5,"Name","conv_1")
    reluLayer("Name","relu_1")
    flattenLayer("Name","flatten")];
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    imageInputLayer([M M 1], "Name", "imageinput_1", "Normalization", "none")
    convolution2dLayer([3 3],10,"Name","conv_2")
    reluLayer("Name","relu_2")
    convolution2dLayer([2 2],5,"Name","conv_3")
    reluLayer("Name","relu_3")
    flattenLayer("Name","flatten_1")];
lgraph = addLayers(lgraph,tempLayers);

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
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Connect Layer Branches

lgraph = connectLayers(lgraph,"flatten","depthcat/in1");
lgraph = connectLayers(lgraph,"flatten_1","depthcat/in2");

%% Plot Layers

figure;
plot(lgraph);