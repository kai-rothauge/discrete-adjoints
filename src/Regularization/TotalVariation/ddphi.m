function y = ddphi(x,eps,varargin)

if nargin > 2
    if strcmpi(varargin{1},'huber')
        y = x;
        y(abs(x) <= eps) = 1;
        y(abs(x) > eps)  = 0;
        return
    elseif strcmpi(varargin{1},'pseudo huber') || strcmpi(varargin{1},'pseudo-huber')
        y = eps*(-x.^2./(eps^2 + x.^2) + 1)./sqrt(eps^2 + x.^2);
        return
    end
end

y = eps*(-x.^2./(eps^2 + x.^2) + 1)./sqrt(eps^2 + x.^2);
