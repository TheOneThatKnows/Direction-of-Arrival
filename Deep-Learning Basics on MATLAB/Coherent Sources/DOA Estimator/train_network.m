%%

numObservations = size(labels, 1);

[idxTrain,idxTest] = trainingPartitions(numObservations,[0.7 0.3]);

XTrain = features(:, :, :, idxTrain);
YTrain = labels(idxTrain, :);
XTest = features(:, :, :, idxTest);
YTest = labels(idxTest, :);

%%

numObservations = size(YTrain, 1);

[idxTrain,idxValidation] = trainingPartitions(numObservations,[0.85 0.15]);

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