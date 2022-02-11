
% Regularization settings

regularization_choice = 'TV';      % T  - Tikhonov regularization
                                   % TV - TV regularization
beta = 0;           % Regularization parameter

if strcmpi(reg_choice, 'T')
    W = [];
elseif strcmpi(reg_choice, 'TV')
    phi = '2';      % 1 - phi(x) = sqrt(x^2 + eps)
                    % 2 - phi(x) = Huber norm
    eps = 0.01;
end