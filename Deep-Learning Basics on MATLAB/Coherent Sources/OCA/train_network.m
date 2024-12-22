%%

numObservations = size(labels, 1);

split_idx = round(numObservations * 0.7);
idxTrain = 1:split_idx;
idxTest = split_idx+1:numObservations;

XTrain = features(:, :, :, idxTrain);
YTrain = labels(idxTrain, :);
XTest = features(:, :, :, idxTest);
YTest = labels(idxTest, :);

%%

numObservations = size(YTrain, 1);

split_idx = round(numObservations * 0.85);
idxTrain = 1:split_idx;
idxValidation = split_idx+1:numObservations;

XValidation = XTrain(:, :, :, idxValidation);
YValidation = YTrain(idxValidation, :);

XTrain = XTrain(:, :, :, idxTrain);
YTrain = YTrain(idxTrain, :);

%%

miniBatchSize  = 128;
validationFrequency = floor(size(YTrain, 1)/miniBatchSize);

options = trainingOptions("sgdm", ...
    MiniBatchSize=miniBatchSize, ...
    InitialLearnRate=1e-3, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.1, ...
    LearnRateDropPeriod=20, ...
    Shuffle="every-epoch", ...
    ValidationData={XValidation, YValidation}, ...
    ValidationFrequency=validationFrequency, ...
    Plots="training-progress", ...
    Metrics="rmse", ...
    Verbose=false);

%%

net = trainnet(XTrain,YTrain,lgraph,"mse",options);
% net = trainnet(XTrain,YTrain,net,"mse",options);

%% Test Network

numObservationsTest = numel(XTest);
for i=1:numObservationsTest
    sequence = XTest{i};
    sequenceLengthsTest(i) = size(sequence,1);
end

[sequenceLengthsTest,idx] = sort(sequenceLengthsTest);
XTest = XTest(idx);
TTest = TTest(idx);

acc = testnet(net,XTest,TTest,"accuracy")

scores = minibatchpredict(net,XTest);
YTest = scores2label(scores,classNames);

figure
confusionchart(TTest,YTest)