classdef MyDataStore < matlab.io.Datastore & matlab.io.datastore.Shuffleable
    properties
        Features1
        Features2
        Labels
    end
    
    properties (Access = private)
        NumObservations
        CurrentIndex
    end
    
    methods
        function ds = MyDataStore(features_1, features_2, labels)
            ds.Features1 = features_1;
            ds.Features2 = features_2;
            ds.Labels = labels;
            ds.NumObservations = size(labels, 2);
            ds.CurrentIndex = 1;
        end
        
        function tf = hasdata(ds)
            tf = ds.CurrentIndex <= ds.NumObservations;
        end
        
        function [data, info] = read(ds)
            if hasdata(ds)
                data{1} = ds.Features1(:, :, :, ds.CurrentIndex);
                data{2} = ds.Features2(:, ds.CurrentIndex);
                data{3} = ds.Labels(:, ds.CurrentIndex);
                info = struct();
                ds.CurrentIndex = ds.CurrentIndex + 1;
            else
                data = [];
                info = struct();
            end
        end

        function ds = shuffle(ds)
            shuffledIndices = randperm(ds.NumObservations);
            ds.Features1 = ds.Features1(:, :, :, shuffledIndices);
            ds.Features2 = ds.Features2(:, shuffledIndices);
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