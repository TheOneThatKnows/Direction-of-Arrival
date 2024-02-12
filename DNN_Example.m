trainImagesFile = "train-images-idx3-ubyte.gz";
testImagesFile = "t10k-images-idx3-ubyte.gz";

XTrain = processImagesMNIST(trainImagesFile);


%%
clear; clc; close all;

% Generate random data for demonstration
X_train = rand(1000, 20);
Y_train = randi([0, 1], 1000, 1);
X_test = rand(200, 20);
Y_test = randi([0, 1], 200, 1);

% Convert responses to categorical format
Y_train = categorical(Y_train);
Y_test = categorical(Y_test);

% Define the model architecture
layers = [
    imageInputLayer([20 1 1])
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(2)
    softmaxLayer    % Use softmax activation for classification
    classificationLayer
];

% Compile the model
options = trainingOptions('adam', ...
    'MaxEpochs', 10, ...
    'MiniBatchSize', 32);

% Train the model
net = trainNetwork(X_train', Y_train, layers, options);

% Evaluate the model
Y_pred = classify(net, X_test');
accuracy = sum(Y_pred == Y_test) / numel(Y_test);
disp(['Test accuracy: ', num2str(accuracy * 100), '%']);
