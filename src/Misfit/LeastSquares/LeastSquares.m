classdef LeastSquares < Misfit

    properties (SetAccess = private)

    end
    
    methods (Access = public)
        function this = LeastSquares(varargin)
            this.name = 'L2';
        end
        
        function m = value(this,d_sim)
            m = this.M(d_sim);
        end 
        
        function m = M(this,d_sim)
            res_d = d_sim - this.d_obs;
            m = 0.5*res_d(:)'*res_d(:);
        end
        
        function dm = dM(this,d_sim)
            res_d = d_sim - this.d_obs;
            dm = res_d;
        end
        
        function ddm = ddM(this,d_sim,varargin)
            if nargin > 2
                ddm = varargin{1};
            else
                ddm = speye(numel(d_sim));
            end
        end
    end
end
