classdef Grid < handle
    
    properties (SetAccess = protected)
        
        dim  = 2;            	% Number of spatial dimensions

        Grad = [];
        Div  = [];
        Curl = [];
        Lap  = [];
        
        % Figure and axis handles for display purposes
        handles = [];
    end
    
    methods (Access = public)
        function this = Grid(varargin)
            
        end

        function set(this,varargin)
            
            if nargin > 1
                for i = 1:2:nargin-1
                    parameter = varargin{i};

                    if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                        disp('  ''dim''');
                    elseif strcmpi(parameter,'dim')
                        this.set_dim(varargin{i+1});
                    end
                end
            end
        end
    end

    methods (Access = private)
        
        % ----------------------- Mutator functions -----------------------

        function flag = set_dim(this,dim)
            
            if dim < 1 || dim > 3
                flag = 1;
            else
                this.dim = dim;
                
                flag = 0;
            end
        end
        
        % -----------------------------------------------------------------

        
    end
end
