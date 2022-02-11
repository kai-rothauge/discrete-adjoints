
% ==================== Minimization procedure settings ====================
%
% --------------------- Minimization procedure codes ----------------------
%
% 'SD':    Steepest descent
% 'NCG':   Nonlinear CG
% 'GN':    Gauss-Newton
% 'LM':    Levenberg-Marquardt
% 'N':     Newton
% 'BFGS':  BFGS
%
% ------------------- Minimization procedure parameters -------------------
%
% Steepest decent:
%
%     line_search_settings (see below)
%
% Nonlinear CG:
% 
%     variant           Nonlinear conjugate gradient variant:
%                       1 - Fletcher Reeves
%                       2 - Polak Ribiere
%                       3 - Hestenes-Stiefel
%     line_search_settings (see below)
%
% Gauss-Newton:
%
%     solver_settings (see below)
%     line_search_settings (see below)
%
% Levenberg-Marquardt:
% 
%     lambda
%     solver_settings (see below)
%     line_search_settings (see below)
%
% Newton:
%
%     solver_settings (see below)
%     line_search_settings (see below)
%         
% BFGS:
%
%     line_search_settings (see below)
%
% LBFGS:
%
%     max_store         Maximum number of stored vectors   
%     scaling           Scaling used in LBFGS
%                       1 - none
%                       2 - tau_1 = s'*H_0*s/s'*y
%                       3 - tau_2 = -alpha*grad'*s/s'*y
%                       4 - tau_3 = s'*grad/grad*(H_0\y)
%                       5 - use sigma = s'*y/y*(H_0\y) if tau_1 > 1
%                       6 - use sigma = s'*y/y*(H_0\y) if tau_2 > 1
%                       7 - use sigma = s'*y/y*(H_0\y) if tau_3 > 1
%     line_search_settings (see below)
%
% -------------------------------------------------------------------------
% ============================ Solver settings ============================
%
% ------------------------------ Solver codes -----------------------------
%
% 'B':    Matlab 'backslash'
% 'CG':   Conjugate gradient
%
% ----------------------------- Solver settings ---------------------------
%
% Backslash:
%
%     code              Solver code
%
% Conjugate gradient:
%     
%     code              Solver code
%     err_tol           Relative error tolerance
%     max_iters         Maximum number of inner iterations
%
%
% -------------------------------------------------------------------------
% ========================= Line search settings ==========================
% 
%     c_1 
%     c_2  
%     alpha_max         Maximum step length
%     max_iters         Maximum number of line search iterations
%
% -------------------------------------------------------------------------
% ========================== Switching criteria ===========================
%
% The following criteria are available for switching minimization
% procedures:
%
%     rel_tol           Switch once the norm of the current gradient has
%                       reached the value of 'rel_tol' relative to the norm 
%                       of the first gradient in the current procedure
%     max_iters         Switch after this many iterations of the current
%                       procedure
%
% The nth minimization procedure is switched to the (n+1)th procedure in the 
% list, where n+1 is the positive integer given in the 'procedure' field.
% The list does not have to be traversed sequentially after the first procedure.
%
% ========================= Minimization settings =========================

% Global minimization settings
rel_tol   = 1.0e-5;         % Relative tolerance
max_iters = 40;             % Maximum number of iterations

% ---------------------- Minimization procedure list ----------------------

% 1st procedure:

    procedure.code = 'SD';            % Steepest descent
    procedure.color = 'b';

    % Line search settings
    procedure.line_search.c_1         = 1.0e-4;
    procedure.line_search.c_2         = 0.45;
    procedure.line_search.alpha_max   = 1e10;
    procedure.line_search.max_iters   = 10;

    % Switching:
        procedure.switch{1}.max_iters = 2;
        procedure.switch{1}.procedure = 2;

        procedure.switch{2}.rel_tol   = 1.0e-1;
        procedure.switch{2}.procedure = 2;

    procedure_list{1} = procedure;
    clear procedure

% 2nd procedure:

    procedure.code = 'NCG';       % Nonlinear CG
    procedure.variant = 1;

    procedure.color = 'r';

    % Line search settings
    procedure.line_search.c_1         = 1.0e-4;
    procedure.line_search.c_2         = 0.45;
    procedure.line_search.alpha_max   = 1e10;
    procedure.line_search.max_iters   = 10;

    % Switching:
        procedure.switch{1}.max_iters = 2;
        procedure.switch{1}.procedure = 3;

        procedure.switch{2}.rel_tol   = 1.0e-4;
        procedure.switch{2}.procedure = 1;

    procedure_list{2} = procedure;
    clear procedure

% 3rd procedure:

    procedure.code = 'GN';       % Gauss-Newton

    procedure.color = 'k';

    % Solver settings
    procedure.solver.code             = 'CG';
    procedure.solver.rel_tol          = 1.0e-3;
    procedure.solver.max_iters        = 15;

    % Line search settings
    procedure.line_search.c_1         = 1.0e-4;
    procedure.line_search.c_2         = 0.45;
    procedure.line_search.alpha_max   = 1e10;
    procedure.line_search.max_iters   = 10;

    % Switching:
        procedure.switch{1}.max_iters = 30;
        procedure.switch{1}.procedure = 1;

        procedure.switch{2}.rel_tol   = 1.0e-2;
        procedure.switch{2}.procedure = 1;

    procedure_list{1} = procedure;
    clear procedure

% =========================================================================

% Save the selected quantities at each minimization iteration
save_history.x                    = 1;      % Solution
save_history.grad_f               = 1;      % Gradient of objective function
save_history.f                    = 1;      % Objective function value
save_history.norm_grad_f          = 1;      % Norm of gradient of objective function value
save_history.alpha                = 1;      % Line search step length
save_history.cpu_time             = 1;      % Total CPU time so far

% Figure numbers for plotting (set to 0 to switch off plotting)
plot_settings.fig_num.x           =  0;    % Current solution
plot_settings.fig_num.grad_f      =  0;    % Current gradient of objective function
plot_settings.fig_num.history     =  99;    % History of selected quantities
plot_settings.fig_num.line_search = 100;    % Line search results

% Plot the history of the selected quantities at each minimization iteration, 
% does not plot quantities that are not saved in history
plot_settings.f                   = 1;      % Objective function value
plot_settings.norm_grad_f         = 1;      % Norm of gradient of objective function value
plot_settings.alpha               = 1;      % Line search step length
plot_settings.cpu_time            = 1;      % Total CPU time so far

% =========================================================================
