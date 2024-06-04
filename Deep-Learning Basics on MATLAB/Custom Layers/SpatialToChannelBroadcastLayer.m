classdef SpatialToChannelBroadcastLayer < nnet.layer.Layer & nnet.layer.Formattable
    methods
        function layer = SpatialToChannelBroadcastLayer(name, opt)
            layer.Name = name;
            layer.Description = "Spatial-to-Channel Broadcast Layer";
        end

        function SSCB = predict(layer, data)
            % data has dimensions [C, B] (channel, batch)
            % Add a spatial dimension to the data
            % SSCB = dlarray(data, "SCB");
            SSCB = dlarray(data, "SSCB");
        end
    end
end
