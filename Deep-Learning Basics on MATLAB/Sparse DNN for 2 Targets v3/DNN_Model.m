%% Create Layer Graph

lgraph = layerGraph();

%% Add Layer Branches

N = 10; Q = 121;

tempLayers = [
    featureInputLayer(2*N-1, "Normalization", "none")
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