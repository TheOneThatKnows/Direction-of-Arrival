%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

M = 5; Q = 121;

tempLayers = [
    featureInputLayer(M*M, "Normalization", "none")
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(150)
    reluLayer
    fullyConnectedLayer(Q)
    regressionLayer];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

%% Plot Layers

figure;
plot(lgraph);