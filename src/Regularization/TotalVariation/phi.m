function y = phi(x,eps,varargin)

x = abs(x);

if nargin > 2
    if strcmpi(varargin{1},'huber')
        y = x;
        y(x <= eps) = 0.5*(x(x <= eps).^2);
        y(x > eps)  = eps*(x(x > eps) - 0.5*eps);
        return
    elseif strcmpi(varargin{1},'pseudo huber') || strcmpi(varargin{1},'pseudo-huber')
        y = eps^2*(sqrt(1 + (x/eps).^2) - 1); 
        return
    end
end

y = eps*sqrt(eps^2 + x.^2); 
