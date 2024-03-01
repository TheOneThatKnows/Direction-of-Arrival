classdef CRN_DataStore < matlab.io.Datastore
    properties
        Features
        Labels
    end
    
    properties (Access = private)
        NumObservations
        CurrentIndex
    end
    
    methods
        function ds = CRN_DataStore(features, labels)
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
                % data = struct('Features', ds.Features(:, :, :, ds.CurrentIndex), ...
                %               'Labels', ds.Labels(:, ds.CurrentIndex));
                data{1} = ds.Features(:, :, :, ds.CurrentIndex);
                data{2} = ds.Labels(:, ds.CurrentIndex);
                info = struct();
                ds.CurrentIndex = ds.CurrentIndex + 1;
            else
                data = [];
                info = struct();
            end
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