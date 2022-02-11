function s = settings_SwiftHohenburg(option)

if strcmpi(option,'PDE')
    
    s.PDE.name = 'SwiftHohenburg';
                 
elseif strcmpi(option,'grid')

    % Domain is periodic
    s.grid.x.max = 40*pi;
    s.grid.x.min = 0;
    s.grid.x.N   = 2^5;
    
    s.grid.y.max = 40*pi;
    s.grid.y.min = 0;
    s.grid.y.N   = 2^5;
    
    % Spatial differentiation order
    s.grid.sd = 6;          % has to be an even integer:
                          % 0 (pseudospectral), 2, 4, 6
                          % Do not select '0' when using a Rosenbrock-type
                          % method

elseif strcmpi(option,'integrator')
    
    s.integrator.choice = 'ETDRK';
    
    s.integrator.scheme = 'Krogstad';   % ETD scheme ('Euler', 'Cox-Matthews', 
                                 %     'Krogstad', 'Hochbruck-Ostermann')
    s.integrator.rosenbrock = 0;      % Use Rosenbrock-type approach (ignored when
                                    % using a pseuospecteal method)
                                    
    s.integrator.max_dt = 0.25;
    
    s.integrator.contour.shape = 'parabolic';   % Contour shape for contour integration,
                                     %      'parabolic' or 'hyperbolic'
    s.integrator.contour.M = 2^5;           % Number of quadrature nodes for contour integral
    
elseif strcmpi(option,'observations')
    
    obs_time_int   = 1;
    obs_time_start = 1;
    obs_time_end   = 20.0;
    
    s.observations.times = obs_time_start:obs_time_int:obs_time_end;
    
    s.observations.noise = 2;         % Noise level as percentage
    
elseif strcmpi(option,'objective') || strcmpi(option,'objective function')
    
    % misfit settings
    s.misfit.choice = 'L2';                         % Least squares misfit
    
    % Regularization settings
    s.regularization{1}.choice  = 'TV';             % Total variation regularization
    s.regularization{1}.beta    = 20;                % Regularization parameter
    s.regularization{1}.norm    = 'pseudo huber';   % 'huber' or 'pseudo huber'
    s.regularization{1}.epsilon = 0.05;              % Parameter needed by Huber norm
    
    s.regularization{2}.choice  = 'Tikhonov';       % Tikhonov regularization
    s.regularization{2}.beta    = 0;              % Regularization parameter
    
elseif strcmpi(option,'minimization')

    % Stopping criteria
    s.rel_tol   = 1e-6;             % Relative tolerance for norm of gradient 
    s.max_iters = 1000;               % Maximum number of minimization iterations
    
    s.procedure.name = 'NCG';        % 'SD':    Steepest descent
                                    % 'NCG':   Nonlinear CG
                                    % 'GN':    Gauss-Newton
                                    % 'BFGS':  L-BFGS
    
    if strcmpi(s.procedure.name,'NCG')  % Nonlinear conjugate gradient settings
        s.procedure.variant = 1;        % NCG variant:
                                        % 1: Fletcher Reeves
                                        % 2: Polak-Ribiere
                                        % 3: Hestenes-Stiefel
    elseif strcmpi(s.procedure.name,'BFGS')     % L-BFGS settings
        s.procedure.M = 20;             % Number of past updates to keep in storage
    end

    % Line search settings
    s.line_search.c_1       = 1.0e-4;   % Parameters for strong Wolfe conditions
    s.line_search.c_2       = 0.45;     %                  
    s.line_search.alpha_max = 1e2;      % Maximum step length
    s.line_search.max_iters = 30;       % Maximum number of line search iterations

elseif strcmpi(option,'display')        % Display/visualization settings

    if strcmpi(varargin{1},'forward')           % Displayed quantities for forward simulation
        s.forward.solution                 = 1;
        s.forward.IC                       = 1;
        s.forward.parameters               = 1;
    elseif strcmpi(varargin{1},'inverse')       % Displayed quantities for parameter estimation
        %   Parameters
        s.inverse.parameters.actual        = 0;  
        s.inverse.parameters.current       = 1;
        s.inverse.parameters.initial       = 0;
        %   Data measurements (final time step only)
        s.inverse.data.actual              = 0;
        s.inverse.data.initial             = 0;
        s.inverse.data.current             = 0;
        s.inverse.data.difference          = 0;
        %   Solution (forward and adjoint), displayed during integration
        s.inverse.solution.forward         = 0;
        s.inverse.solution.adjoint         = 0;
        %   Minimization data
        s.inverse.minimization.history     = 1;
        s.inverse.minimization.line_search = 1;
    end
    
elseif strcmpi(option,'diagnostics')         % Diagnostics settings
    
    % Tests for integration scheme
    s.integration.derivative.DtDy     = 1;
    s.integration.derivative.DtDy_inv = 1;
    s.integration.derivative.DtDm     = 1;
    s.integration.derivative.DyDm     = 1;
    
    s.integration.transpose.DtDy      = 1;      
    s.integration.transpose.DtDy_inv  = 1;   % This is the adjoint test
    s.integration.transpose.DtDm      = 1; 
    s.integration.transpose.DyDm      = 1; 

    s.integration.inverse.t           = 1;
    s.integration.inverse.DtDy        = 1;
    s.integration.inverse.DtDyT       = 1;
    
    % Tests for objective function
    s.objective.gradient          = 0;
    
    % Tests for misfit
    s.misfit.gradient             = 0;
    s.misfit.derivative           = 0;
    
    % Test for regularization
    s.regularization.derivative   = 0;

end