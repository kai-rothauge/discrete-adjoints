classdef Grid_Point < handle
    
    properties (SetAccess = public)
        
        dim = 2;
        x = [];     % Coordinates of grid point
    end
    
    methods (Access = public)
        function this = Grid_Axis(varargin)
            
            this.dim = nargin;
            for i = 1:nargin
                x(i) = varargin{i};
            end
        end
    end
end
