function y = dphi(x,eps,varargin)

if nargin > 2
    if strcmpi(varargin{1},'huber')
        y = x;
        y(abs(x) <= eps) = x(abs(x) <= eps);
        y(abs(x) > eps)  = eps*(sign(x(abs(x) > eps)));
        return
    elseif strcmpi(varargin{1},'pseudo huber') || strcmpi(varargin{1},'pseudo-huber')
        y = eps*x./sqrt(eps^2 + x.^2);
        return
    end
end

y = eps*x./sqrt(eps^2 + x.^2);
