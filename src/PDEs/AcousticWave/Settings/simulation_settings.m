
% Simulation settings

% Select partial differential equation
PDE = 'AWE';

% ============================ Simulation mode ============================

problem = 1;             % 1: Forward problem
                         % 2: Inverse problem
                         % 3: Diagnostics

% ============================= Settings files ============================
     
% Names of settings files
file_names.imaging            = 'Input/Settings/imaging_settings.m';
file_names.solver.PDE         = 'Input/Settings/PDE_settings.m';
file_names.solver.integration = 'Input/Settings/integration_settings.m';

if problem == 2 || problem == 3
    file_names.misfit         = 'Input/Settings/misfit_settings.m';
    file_names.minimization   = 'Input/Settings/minimization_settings.m';
end

if problem == 3
    file_names.diagnostic     = 'Input/Settings/diagnostic_settings.m';
end

% ============================ Output settings ============================

verbose = 1;
padding = '    ';

