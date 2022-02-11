classdef Tikhonov < Regularization

    properties (SetAccess = private)

        W     = [];        
        m_ref = [];
    end
    
    methods (Access = public)
        function this = Tikhonov(varargin)
            this.name = 'Tikhonov';

            if nargin == 1 && isstruct(varargin{1})
                this.beta = varargin{1}.beta;    % Regularization parameter
            end
            
            if nargin > 1
                this.set(varargin{:});
            end
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this,varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter, 'm_ref')
                    this.set_m_ref(varargin{i+1});
                elseif strcmpi(parameter,'W')
                    this.set_W(varargin{i+1});
                end
            end
        end
        
        function set_W(this,W)
            this.W = W;
        end
        
        function set_m_ref(this,m_ref)
            this.m_ref = m_ref;
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function out = get(this,parameter)
            out = [];
            
            if strcmpi(parameter,'m_ref')
                out = this.m_ref;
            elseif strcmpi(parameter,'W')
                out = this.W;
            end
        end
        
        function out = get_m_ref(this)
            out = this.m_ref;
        end
        
        function out = get_W(this)
            out = this.W;
        end
        
        % -----------------------------------------------------------------
        
        function r = value(m,varargin)
            r = this.R(m);
        end
        
        function r = R(this,m)
            n_m_ref = 1;
            
            if isempty(this.W)
                res_m = m(:) - this.m_ref(:);
            else
                res_m = this.W*(m(:) - this.m_ref(:));
            end
            r = 0.5*n_m_ref*(res_m'*res_m);
            r = this.beta*r;
        end
        
        function dr = dR(this,m)
            n_m_ref = 1;
            
            res_m = m(:) - this.m_ref(:);
            if isempty(this.W)
                dr = res_m;
            else
                dr = this.W'*(this.W*res_m);
            end
            dr = n_m_ref*dr;
            dr = this.beta*dr;
        end
        
        function ddr = ddR(this,m,varargin)
            n_m_ref = 1;
            
            if nargin > 2 && ~isempty(varargin{1})
                u = varargin{1};
                if isempty(this.W)
                    ddr = u;
                else
                    ddr = this.W'*(this.W*u);
                end
            else
                if isempty(this.W)
                    ddr = speye(length(m(:)));
                else
                    ddr = this.W'*this.W;
                end
            end
            ddr = n_m_ref*ddr;
            ddr = this.beta*ddr;
        end
        
        function result = run_tests(this)
            
            if isempty(this.m_ref)
                N = 1000;
                this.m_ref = 100*rand(N,1);
            else
                N = numel(this.m_ref);
            end
            result = run_tests@Regularization(this,N);
        end
    end
end
