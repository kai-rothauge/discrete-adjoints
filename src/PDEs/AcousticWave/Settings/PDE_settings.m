
% PDE settings

% ---------------------------- Settings files -----------------------------
     
file_names.domain       = 'Input/Settings/domain_settings.m';

file_names.bulk_modulus = 'Input/Data/kappa_actual.mat';
file_names.density      = 'Input/Data/rho_actual.mat';
% 
% file_names.bulk_modulus = 'Input/Data/kappa_reference.mat';
% file_names.density      = 'Input/Data/rho_reference.mat';

% --------------------- Spatial differentiation order ---------------------

% Spatial differentiation order, has to be an even integer
spatial_differentiation_order = 2;