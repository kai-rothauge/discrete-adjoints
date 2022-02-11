classdef Misfit < handle

    properties (SetAccess = protected)

        name  = '';          % Name of misfit function
        sim   = [];          % Link to Simulation object
        d_obs = [];          % Observed data
    end
    
    methods (Access = public)
        function this = Misfit(varargin)
            
        end
        
        function set(this,varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter, 'd_obs')
                    this.set_d_obs(varargin{i+1});
                elseif strcmpi(parameter,'sim')
                    this.sim = varargin{i+1};
                end
            end
        end
        
        function set_d_obs(this,d_obs)
            this.d_obs = d_obs;
        end
        
        % -----------------------------------------------------------------
        
        function out = get(this,parameter)
            
            if strcmpi(parameter, 'd_obs')
                out = this.d_obs;
            elseif strcmpi(parameter,'sim')
                out = this.sim;
            end
        end
        
        % -----------------------------------------------------------------
        
        function grad_m = grad_M(this,d_sim)
            grad_m = this.gradient(d_sim);
        end
        
        function grad_m = gradient(this,d_sim)
            grad_m = this.sim.J_times(this.dM(d_sim),'transpose',1);
        end
        
        % -------------------------- Diagnostics --------------------------
        
        function result = run_tests(this)
            
            test_gradient = 1;
            
            this.sim.generate_observations();
            this.sim.PDE.set('parameters','initial guess');
            
            result = [];

            if test_gradient
                m  = this.sim.get('m');
                dm = this.sim.get('dm');
                
                d = this.sim.generate_data(m);
                
                M      = this.M(d);
                grad_M = this.grad_M(d);
                
                for i = 1:10
                    h = 1*0.1^(i+1);

                    m_pert = this.M(this.sim.generate_data(m + h*dm));

                    err_M(i)      = m_pert - M;
                    err_grad_M(i) = m_pert - M - h*grad_M'*dm;
                end
                this.sim.PDE.set('m',m);
                
                err_M = err_M(:); err_M
                err_grad_M = err_grad_M(:); err_grad_M
            end
        end
    end
    
    methods (Access = public, Abstract)
        
        m   = value(this,d_sim);
        m   = M(this,d_sim);
        dm  = dM(this,d_sim);
        ddm = ddM(this,d_sim,varargin);
    end
end
