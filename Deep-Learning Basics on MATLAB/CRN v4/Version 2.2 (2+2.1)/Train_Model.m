%%

split_percent = 0.7;
split_idx = size(labels, 2) * split_percent;

X_Train = features(:, :, :, 1:split_idx);
Y_Train = labels(:, 1:split_idx);
% X_Test = features(:, :, :, split_idx+1:end);
% Y_Test = labels(:, split_idx+1:end);
X_Test = features(:, :, :, split_idx+1:idx);
Y_Test = labels(:, split_idx+1:idx);

dsTrain = MyDataStore(X_Train, Y_Train);
dsTest = MyDataStore(X_Test, Y_Test);

%%

% Compile the model with specified loss function
options = trainingOptions('adam', ...
    'InitialLearnRate', 0.001, ...
    'MaxEpochs', 12, ...
    'MiniBatchSize', 128, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress');

%%

net = trainNetwork(dsTrain, lgraph, options);