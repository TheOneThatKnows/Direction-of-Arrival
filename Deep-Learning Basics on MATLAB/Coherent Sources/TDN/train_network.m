%%

numObservations = size(labels, 2);

[idxTrain,idxTest] = trainingPartitions(numObservations,[0.7 0.3]);

XTrain = features(:, :, :, idxTrain);
YTrain = labels(idxTrain);
XTest = features(:, :, :, idxTest);
YTest = labels(idxTest);

%%

numObservations = size(YTrain, 2);

[idxTrain,idxValidation] = trainingPartitions(numObservations,[0.85 0.15]);

XValidation = XTrain(:, :, :, idxValidation);
YValidation = YTrain(idxValidation);

XTrain = XTrain(:, :, :, idxTrain);
YTrain = YTrain(idxTrain);

%%

dsTrain = MyDataStore(XTrain, YTrain);
dsValid = MyDataStore(XValidation, YValidation);
dsTest = MyDataStore(XTest, YTest);

%%

options = trainingOptions("sgdm", ...
    InitialLearnRate=0.001, ...
    MaxEpochs=2, ...
    Shuffle="every-epoch", ...
    ValidationData=dsValid, ...
    ValidationFrequency=30, ...
    Plots="training-progress", ...
    Metrics="accuracy", ...
    Verbose=false);

%%

% net = trainnet(dsTrain,lgraph,"crossentropy",options);
net = trainnet(dsTrain,net,"crossentropy",options);

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