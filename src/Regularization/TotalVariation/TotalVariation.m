classdef TotalVariation < Regularization

    properties (SetAccess = private)
        
        norm    = 'pseudo huber';
        epsilon = 0.1;
        grad    = [];
    end
    
    methods (Access = public)
        function this = TotalVariation(varargin)
            this.name = 'TV';

            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
                this.beta    = s.beta;          % Regularization parameter
                this.norm    = s.norm;          % 'huber' or 'pseudo huber'
                this.epsilon = s.epsilon;       % Parameter needed by Huber norm
            end
            
            if nargin > 1
                this.set(varargin{:});
            end
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this,varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'norm')
                    this.set_norm(varargin{i+1});
                elseif strcmpi(parameter,'eps') || strcmpi(parameter,'epsilon')
                    this.set_epsilon(varargin{i+1});
                elseif strcmpi(parameter,'grad') || strcmpi(parameter,'G')
                    this.set_grad(varargin{i+1});
                end
            end
        end
        
        function set_norm(this,norm)
            this.norm = norm;
        end
        
        function set_epsilon(this,epsilon)
            this.epsilon = epsilon;
        end
        
        function set_grad(this,grad)
            this.grad = grad;
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function out = get(this,parameter)
            out = [];
            
            if strcmpi(parameter,'norm')
                out = this.norm;
            elseif strcmpi(parameter,'eps') || strcmpi(parameter,'epsilon')
                out = this.epsilon;
            elseif strcmpi(parameter,'grad') || strcmpi(parameter,'G')
                out = this.grad;
            end
        end
        
        function out = get_norm(this)
            out = this.norm;
        end
        
        function out = get_epsilon(this)
            out = this.epsilon;
        end
        
        function out = get_grad(this)
            out = this.grad;
        end
        
        % -----------------------------------------------------------------

        function r = value(m,varargin)
            r = this.R(m);
        end
        
        function r = R(this,m)
            G = this.grad;
            r = sum(this.phi(G*m(:)));
            r = this.beta*r;
        end
        
        function dr = dR(this,m)
            G  = this.grad;
            dr = G'*this.dphi(G*m(:));
            dr = this.beta*dr;
        end
        
        function ddr = ddR(this,m,varargin)
            G = this.grad;
            if nargin > 2 && ~isempty(varargin{1})
                u = varargin{1};
                ddr = G'*(this.ddphi(G*m(:)).*(G*u(:)));
            else
                Gm = this.ddphi(G*m(:));
                ddr = G'*spdiags(Gm,0,length(Gm),length(Gm))*G;
            end
            ddr = this.beta*ddr;
        end
        
        function result = run_tests(this)
            
            N = size(this.grad,2);
            
            x = linspace(-1,1,N+1); x = x(:);
            
%             eps = 0.1;
%             
%             figure(11)
%             
%             ax = abs(x); y1 = ax;
%             y1(ax <= eps) = (0.5/eps)*(ax(ax <= eps).^2);
%             y1(ax > eps)  = ax(ax > eps) - 0.5*eps;   
%             y2 = sqrt(eps^2 + ax.^2) - eps;
%             subplot(2,2,1), plot(y1), hold on, plot(y2), hold off
% 
%             dy1 = x;
%             dy1(abs(x) <= eps) = (1/eps)*x(abs(x) <= eps);
%             dy1(abs(x) > eps)  = (sign(x(abs(x) > eps)));
%             dy2 = x./sqrt(eps^2 + x.^2);
%             subplot(2,2,2), plot(dy1), hold on, plot(dy2), hold off
% 
%             ddy1 = x;
%             ddy1(abs(x) <= eps) = 1/eps;
%             ddy1(abs(x) > eps)  = 0;
%             ddy2 = (-x.^2./(eps^2 + x.^2) + 1)./sqrt(eps^2 + x.^2);
%             subplot(2,2,3), plot(ddy1), hold on, plot(ddy2), hold off
            
            result = run_tests@Regularization(this,N);
        end
    end
    
    methods (Access = private)
        function y = phi(this,x)

            eps = this.epsilon;
            x   = abs(x);
            
            if strcmpi(this.norm,'huber')
                y = x;
                y(x <= eps) = (0.5/eps)*(x(x <= eps).^2);
                y(x > eps)  = (x(x > eps) - 0.5*eps);
            elseif strcmpi(this.norm,'pseudo huber') || strcmpi(this.norm,'pseudo-huber')
                y = sqrt(eps^2 + x.^2) - eps;
            else
                y = sqrt(eps^2 + x.^2);
            end
        end
        
        function y = dphi(this,x)

            eps = this.epsilon;

            if strcmpi(this.norm,'huber')
                y = x;
                y(abs(x) <= eps) = (1/eps)*x(abs(x) <= eps);
                y(abs(x) > eps)  = (sign(x(abs(x) > eps)));
            else
                y = x./sqrt(eps^2 + x.^2);
            end
        end
        
        function y = ddphi(this,x)

            eps = this.epsilon;

            if strcmpi(this.norm,'huber')
                y = x;
                y(abs(x) <= eps) = 1/eps;
                y(abs(x) > eps)  = 0;
            else
                y = (-x.^2./(eps^2 + x.^2) + 1)./sqrt(eps^2 + x.^2);
            end
        end
    end
end
