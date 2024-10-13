%%

split_percent = 0.7;
split_idx = size(labels, 2) * split_percent;

X_Train_1 = features_1(:, :, :, 1:split_idx);
X_Train_2 = features_2(:, 1:split_idx);
Y_Train = labels(:, 1:split_idx);
% X_Test = features(:, :, :, split_idx+1:end);
% Y_Test = labels(:, split_idx+1:end);
X_Test_1 = features_1(:, :, :, split_idx+1:idx);
X_Test_2 = features_2(:, split_idx+1:idx);
Y_Test = labels(:, split_idx+1:idx);

dsTrain = MyDataStore(X_Train_1, X_Train_2, Y_Train);
dsTest = MyDataStore(X_Test_1, X_Test_2, Y_Test);

%% Opening up space in RAM

clear  features_1 features_2 lables X_Train_1 X_Train2 Y_Train X_Test_1 X_Test_2 Y_Test dsTest

%%

% Compile the model with specified loss function
options = trainingOptions('adam', ...
    'InitialLearnRate', 0.0005, ...
    'MaxEpochs', 12, ...
    'MiniBatchSize', 128, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress');
% 'ValidationData', dsTest

%%

net = trainNetwork(dsTrain, lgraph, options);