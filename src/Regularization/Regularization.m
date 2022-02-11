classdef Regularization < handle

    properties (SetAccess = public)

        name = '';          % Name of regularization method
        beta = 0;           % Regularization parameter
    end
    
    methods (Abstract)
        r   = value(this,m,varargin);
        r   = R(this,m,varargin);
        dr  = dR(this,m,varargin);
        ddr = ddR(this,m,varargin);
    end
    
    methods (Access = public)
        function this = Regularization()

        end
        
        % -------------------------- Diagnostics --------------------------
        
        function result = run_tests(this,N)

            result = [];
            
            disp(['Testing ' this.name ' Regularization:'])

            m  = 100*rand(N,1);
            dm = 10*rand(N,1);
            
            r   = this.R(m);
            dr  = this.dR(m);
            ddr = this.ddR(m,dm);

            err_R   = zeros(10,1); 
            err_dR  = zeros(10,1); 
            err_ddR = zeros(10,1); 
            
            for i = 1:10
                h = 0.1^(i-2);
                
                r_pert = this.R(m + h*dm);

                err_R(i)   = (r_pert - r)/r;
                err_dR(i)  = (r_pert - r - h*dr'*dm(:))/r;
                err_ddR(i) = (r_pert - r - h*dr'*dm(:) - 0.5*h^2*ddr(:)'*dm(:))/r;
            end
            err_R, err_dR, err_ddR
        end
    end
end
