classdef SSSpatialToChannelBroadcastLayer < nnet.layer.Layer & nnet.layer.Formattable
    methods
        function layer = SSSpatialToChannelBroadcastLayer(name, opt)
            layer.Name = name;
            layer.Description = "Spatial-to-Channel Broadcast Layer";
        end

        function SSCB = predict(layer, data)
            % data has dimensions [C, B] (channel, batch)
            % Add a spatial dimension to the data
            % SSCB = dlarray(data, "SSSCB");
            SSCB = dlarray(data, "CB");
        end
    end
end
