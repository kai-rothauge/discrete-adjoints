classdef Parameter < handle
    
    properties (SetAccess = protected)
        
        name = '';          % Name of parameter
        
        values = [];        % Parameter values
        max = [];           % Maximum value parameter can attain
        min = [];           % Minimum value parameter can attain
        
        N;                  % Number of parameters
        
        % Figure and axis handles for display purposes
        handles = [];
    end
    
    methods (Access = public)
        function this = Parameter(varargin)
            u0 = [];
            qq = [];
        end
    end

    methods (Access = protected)
        
        function D2 = derivative_1D(this,order,N,dx)
            e = ones(N,1);
            if order == 4
                D2 = spdiags([-e 16*e -30*e 16*e -1*e],-2:2,N,N);
                D2(1,end-1:end) = [-1 16]; D2(2,end)   = -1;
                D2(end,1:2)     = [16 -1]; D2(end-1,1) = -1;
                D2 = D2/(12*dx^2);
            elseif order == 6
                D2 = spdiags([2*e -27*e 270*e -490*e 270*e -27*e 2*e],-3:3,N,N);
                D2(1,end-2:end) = [2 -27 270]; D2(2,end-1:end) = [2 -27]; D2(3,end) = 2;
                D2(end,1:3)     = [270 -27 2]; D2(end-1,1:2)   = [-27 2]; D2(end-2,1) = 2;
                D2 = D2/(180*dx^2);
            else            % 2nd order by default
                D2 = spdiags([e -2*e e],-1:1,N,N);
                D2(1,end) = 1; D2(end,1) = 1;
                D2 = D2/dx^2;
            end
        end

    end
end
