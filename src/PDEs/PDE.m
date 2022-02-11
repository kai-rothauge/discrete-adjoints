classdef PDE < handle
    
    properties (SetAccess = protected)
        
        name = '';          % Name of PDE
        
        grid = [];          % Grid object
        m = [];             % Parameter objects
        
        linearity;          % Linearity of PDE (linear, semilinear, or nonlinear)
        
        u0 = [];                 % Initial condition
        q  = [];                 % Source term
        
        % Figure and axis handles for display purposes
        handles = [];
    end
    
    methods (Access = public)
        function this = PDE(varargin)
            
        end
    end
end
