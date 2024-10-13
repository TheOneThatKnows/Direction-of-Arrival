classdef MyDataStore < matlab.io.Datastore & matlab.io.datastore.Shuffleable
    properties
        Features
        Labels
    end
    
    properties (Access = private)
        NumObservations
        CurrentIndex
    end
    
    methods
        function ds = MyDataStore(features, labels)
            ds.Features = features;
            ds.Labels = labels;
            ds.NumObservations = size(labels, 2);
            ds.CurrentIndex = 1;
        end
        
        function tf = hasdata(ds)
            tf = ds.CurrentIndex <= ds.NumObservations;
        end
        
        function [data, info] = read(ds)
            if hasdata(ds)
                data{1} = ds.Features(:, :, 1, ds.CurrentIndex);
                data{2} = ds.Features(:, :, 2, ds.CurrentIndex);
                data{3} = ds.Features(:, :, 3, ds.CurrentIndex);
                data{4} = ds.Labels(:, ds.CurrentIndex);
                info = struct();
                ds.CurrentIndex = ds.CurrentIndex + 1;
            else
                data = [];
                info = struct();
            end
        end

        function ds = shuffle(ds)
            shuffledIndices = randperm(ds.NumObservations);
            ds.Features = ds.Features(:, :, :, shuffledIndices);
            ds.Labels = ds.Labels(:, shuffledIndices);
        end

        function reset(ds)
            ds.CurrentIndex = 1;
        end
        
        function n = get.NumObservations(ds)
            n = ds.NumObservations;
        end

        function output = GetNumberOfObservations(ds)
            output = ds.NumObservations;
        end

        function output = GetCurrentIndex(ds)
            output = ds.CurrentIndex;
        end
    end
end