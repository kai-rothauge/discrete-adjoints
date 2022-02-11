classdef Minimization < handle
    
    properties (SetAccess = public)

        rel_tol = 1e-6;
        max_iters = 10;
        
        cutoff_tol = 0;
        
        procedure = [];

        sim            = [];        % Link to Simulation object
        misfit         = [];        % Link to Misfit object
        regularization = [];        % Link to Regularization object(s)
        line_search    = [];        % Link to Line Search object
        
        % Workspace storage
        m_current = [];
        d_current = [];

        % History of relevant quantities
        M_current = 0;          % Current value of the misfit function
        R_current = 0;          % Current value of the regularization function
        
        grad_M_current = [];    % Current gradient of the misfit function
        grad_R_current = [];    % Current gradient of the regularization function
        
        history = [];
        
        % Display settings
        display = [];
    end
    
    methods (Access = public)
        function this = Minimization(varargin)
            
            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
            else
                s = settings('minimization');
            end
            
            if nargin > 1
                this.set(varargin{:});
            end
    
            % Minimization settings
            this.procedure   = s.procedure;
            this.rel_tol     = s.rel_tol;         % Relative tolerance stopping criterion
            this.max_iters   = s.max_iters;
            if isempty(this.line_search)
                this.line_search = LineSearch(s.line_search);
            end
            this.line_search.set('function',@(m)(this.compute_chi(m)), ...
                                 'gradient',@(m)(this.compute_grad_chi(m)));
            
            % Objective function settings
            s = settings('objective function');
            
            if isempty(this.misfit)
                this.misfit = LeastSquares('sim',this.sim);
            end
            
            if isempty(this.regularization)
                for ri = 1:length(s.regularization)
                    if strcmpi(s.regularization{ri}.choice,'Tikhonov')
                        this.regularization{ri} = Tikhonov(s.regularization{ri});
                    elseif strcmpi(s.regularization{ri}.choice,'TV')
                        this.regularization{ri} = TotalVariation(s.regularization{ri});
                        if ~isempty(this.sim)
                            G = this.sim.get_parameter_gradient();
                            this.regularization{ri}.set('grad',G);
                        end
                    end
                end
            end
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this,varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('Options for ''Minimization::set'':');
                    disp('  ''sim''');
                    disp('  ''rel tol'' or ''relative tolerance''');
                    disp('  ''max iters'' or ''max it'' or ''max iterations''');
                    disp('  ''line search'' or ''ls''');
                    disp('  ''m_ref''');
                    disp('  ''d_obs''');
                    disp('  ''display history''');
                    disp('  ''display line search''');
                elseif strcmpi(parameter,'sim')
                    this.set_sim(varargin{i+1});
                elseif strcmpi(parameter,'rel tol') || strcmpi(parameter,'relative tolerance')
                    this.rel_tol = varargin{i+1};
                elseif strcmpi(parameter,'max iters') || strcmpi(parameter,'max it') || ...
                       strcmpi(parameter,'max iterations')
                    this.max_iters = max(2,varargin{i+1});
                elseif strcmpi(parameter,'line search') || strcmpi(parameter,'ls')
                    this.line_search = varargin{i+1};
                elseif strcmpi(parameter,'d_obs')
                    this.set_d_obs(varargin{i+1});
                elseif strcmpi(parameter,'m_ref')
                    this.set_m_ref(varargin{i+1});
                elseif strcmpi(parameter,'display history')
                    this.display.history.on = varargin{i+1};
                elseif strcmpi(parameter,'display line search')
                    this.display.line_search.on = varargin{i+1};
                end
            end
        end
        
        function set_sim(this,s)
            this.sim = s;
            for ri = 1:length(this.regularization)
                if strcmpi(this.regularization{ri}.name,'TV')
                    G = this.sim.get_parameter_gradient();
                    this.regularization{ri}.set('grad',G);
                end
            end
            this.misfit.set('sim',s);
        end
        
        function set_m_ref(this,m_ref)
            for ri = 1:length(this.regularization)
                this.regularization{ri}.set('m_ref',m_ref);
            end
        end
        
        function set_d_obs(this,d_obs)
            this.misfit.set('d_obs',d_obs);
        end
        
        % ---------------------- Accessor functions -----------------------
        
        function out = get(this,parameter)
            out = [];
            
            if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                disp('Options for ''Minimization::get'':');
                disp('  ''sim''');
                disp('  ''rel tol'' or ''relative tolerance''');
                disp('  ''max iters'' or ''max it'' or ''max iterations''');
                disp('  ''line search'' or ''ls''');
                disp('  ''axes handle''');
            elseif strcmpi(parameter,'sim')
                out = this.sim;
            elseif strcmpi(parameter,'rel tol') || strcmpi(parameter,'relative tolerance')
                out = this.rel_tol;
            elseif strcmpi(parameter,'max iters') || strcmpi(parameter,'max it') || ...
                   strcmpi(parameter,'max iterations')
                out = this.max_iters;
            elseif strcmpi(parameter,'line search') || strcmpi(parameter,'ls')
                out = this.line_search;
            elseif strcmpi(parameter,'axes handle')
                out = this.axes_handle;
            end
        end

        % -----------------------------------------------------------------
        
        function history = run(this,varargin)
            
            this.setup_figures();
            figure(this.display.history.fig)
            this.initialize_history();
            
            if nargin > 1
                m = varargin{1};
            else
                m = this.sim.get('parameters');
            end
            fprintf('Computing current objective function values ........... ');
            chi      = this.compute_chi(m);
            grad_chi = this.compute_grad_chi(m);
            fprintf('done\n');
            
            m_prev        = m;
            grad_chi_prev = grad_chi;
            dm            = [];
            
            norm_grad_chi_0 = norm(grad_chi);
            this.cutoff_tol = this.rel_tol*norm_grad_chi_0;
            
            fprintf('Current values:\n')
            fprintf('    M   = %3.6f     ||grad_M||   = %3.6f\n',this.M_current,norm(this.grad_M_current))
            fprintf('    R   = %3.6f     ||grad_R||   = %3.6f\n',this.R_current,norm(this.grad_R_current))
            fprintf('    chi = %3.6f     ||grad_chi|| = %3.6f\n',chi,norm_grad_chi_0)
            
            this.history.chi(1)      = chi;
            this.history.grad_chi{1} = grad_chi;
            this.history.m{1}        = m;
          	
            for k = 1:this.max_iters
                fprintf('\nIteration %d of %d:\n',k,this.max_iters);
                
                fprintf('Computing search direction ............................ ');
                dm = this.get_search_direction(grad_chi,dm,grad_chi_prev,m-m_prev,k);
                fprintf('done\n');
                
%                 if k > 1 && -grad_chi(:)'*dm(:) <= sqrt(eps)
%                     fprintf('Running tests ......................................... ');
%                     results = this.run_tests(m,dm,chi,grad_chi);
%                     fprintf('done');
%                     keyboard
%                 end

                grad_chi_prev = grad_chi;
                m_prev        = m;
                
                if k > 5 && abs((this.history.chi(k-5)-this.history.chi(k))/this.history.chi(k)) < sqrt(eps)
                    m = m + dm;
                    chi      = this.compute_chi(m);
                    grad_chi = this.compute_grad_chi(m);
                else
                    fprintf('Performing line search ................................ ');
                    t_start = cputime;
                    if k < 4
                        alpha = abs(m(:)'*dm(:)/(dm(:)'*dm(:)));
                    else
                        alpha = 1;
                    end
                    [m,alpha,chi,grad_chi] = this.line_search.run(m,dm,alpha,chi,grad_chi);
                    fprintf('done (%6.2fs, alpha = %6.4e)\n',cputime - t_start,alpha);
                end
                
                m = this.sim.constrain_parameters(m);
                this.sim.set('m',m);
                this.sim.update_plots('m',m);

                norm_grad_chi = norm(grad_chi);
                
                fprintf('Current values:\n')
                fprintf('    M   = %8.6f     ||grad_M||   = %8.6f\n',this.M_current,norm(this.grad_M_current))
                fprintf('    R   = %8.6f     ||grad_R||   = %8.6f\n',this.R_current,norm(this.grad_R_current))
                fprintf('    chi = %8.6f     ||grad_chi|| = %8.6f\n',chi,norm_grad_chi)
                
                % Save and plot history
                this.save_history(k,m,dm);
                if this.display.history.on, this.plot_history; end

                if norm_grad_chi < this.rel_tol*norm_grad_chi_0
                    break;
                end
            end
            
            history = this.history;
        end

        % -------------------------- Diagnostics --------------------------
        
        function result = run_tests(this,varargin)

            result = [];

            m        = varargin{1};
            dm       = varargin{2};
            chi      = varargin{3};
            grad_chi = varargin{4};

            for i = 1:10
                h = 1*0.1^(i+1);

                chi_pert = this.compute_chi(m + h*dm);

                err_chi(i)      = (chi_pert - chi)/chi;
                err_grad_chi(i) = (chi_pert - chi - h*grad_chi'*dm)/chi;
            end

            err_chi      = err_chi(:); err_chi
            err_grad_chi = err_grad_chi(:); err_grad_chi
            
            chi_pert = this.compute_chi(m);
        end
    end

    methods (Access = private)
        
        function dm = get_search_direction(this,grad_chi,varargin)
            
            proc_name = this.procedure.name;
            
            if strcmpi(proc_name,'steepest descent') || strcmpi(proc_name,'SD')
                dm = -grad_chi;
            elseif strcmpi(proc_name,'NCG')
                dm            = varargin{1};
                grad_chi_prev = varargin{2};

                if isempty(dm)
                    dm = -grad_chi;            % Steepest descent at first iteration
                else
                    switch this.procedure.variant
                        case 1          % Fletcher-Reeves
                            beta = (grad_chi'*grad_chi)/(grad_chi_prev'*grad_chi_prev);
                        case 2          % Polak-Ribiere
                            delta_grad_M = grad_chi - grad_chi_prev;
                            beta = max((grad_chi'*delta_grad_M)/(grad_chi_prev'*grad_chi_prev), 0);
                        case 3          % Hestenes-Stiefel
                            delta_grad_M = grad_chi - grad_chi_prev;
                            beta = -(grad_chi'*delta_grad_M)/(delta_grad_M'*dm);
                        case 4          % Dai-Yuan
                            delta_grad_M = grad_chi - grad_chi_prev;
                            beta = -(grad_chi'*grad_chi)/(delta_grad_M'*dm);
                    end

                    dm = -grad_chi + beta*dm;
                end
            elseif strcmpi(proc_name,'Gauss-Newton') || strcmpi(proc_name,'GN')
                m = varargin{3};
                
                GN_fcn = @(x)(this.sim.J_times(this.sim.J_times(x),1) + this.regularization{1}.ddR(m,x));
                
                [dm,flag,relres,iter] = minres(GN_fcn,-grad_chi,1e-6,200);
                iter
                relres
            elseif strcmpi(proc_name,'BFGS') || strcmpi(proc_name,'LBFGS')
                dm            = varargin{1};
                grad_chi_prev = varargin{2};
                s             = varargin{3};
                k             = varargin{4};
                
                M = this.procedure.M;
                
                if ~isfield(this.procedure,'s') || ~isfield(this.procedure,'y') || ...
                        ~isfield(this.procedure,'rho')
                    this.procedure.s   = zeros(length(s),this.max_iters);
                    this.procedure.y   = zeros(length(s),this.max_iters);
                    this.procedure.rho = zeros(this.max_iters,1);
                end

                if isempty(dm)
                    dm = -grad_chi;            % Steepest descent at first iteration
                else
                    this.procedure.s(:,k-1) = s;
                    this.procedure.y(:,k-1) = grad_chi - grad_chi_prev;
                    this.procedure.rho(k-1) = 1/(this.procedure.y(:,k-1)'*this.procedure.s(:,k-1));
                    
                    alpha = zeros(this.max_iters,1);
                    q     = grad_chi;
                    for i = k-1:-1:max(1,k-M)
                        alpha(i) = this.procedure.rho(i)*(this.procedure.s(:,i)'*q);
                        q = q - alpha(i)*this.procedure.y(:,i);
                    end
                
                    y_prev = this.procedure.y(:,k-1);
                    dm = ((this.procedure.s(:,k-1)'*y_prev)/(y_prev'*y_prev))*q;

                    for i = max(1,k-M):k-1
                        beta = this.procedure.rho(i)*(this.procedure.y(:,i)'*dm);
                        dm = dm + (alpha(i) - beta)*this.procedure.s(:,i);
                    end
                    dm = -dm;
                end
            end
            
            % If search direction is not a descent direction because the
            % previous line search failed, switch the sign
            if grad_chi(:)'*dm(:) > -sqrt(eps)
                dm = -grad_chi;
            end
        end
        
        % -----------------------------------------------------------------
        
        function chi = compute_chi(this,m)
            
            m = this.sim.constrain_parameters(m);
            
            if isempty(this.d_current) || isempty(this.m_current) || ...
                                          (norm(m(:)-this.m_current(:)) > 1e-13)
                this.d_current = this.sim.generate_data(m);
                this.sim.display_data(this.d_current);
                this.m_current = m;
            end
            
            % Misfit
            M = this.misfit.M(this.d_current);
            M = M/numel(this.d_current);
            
            % Regularization
            R = 0;
            for ri = 1:length(this.regularization)
                R = R + this.regularization{ri}.R(this.m_current);
            end
            R = R/numel(this.m_current);
            
            this.M_current = M;
            this.R_current = R;
            
            chi = M + R;
        end
        
        function grad_chi = compute_grad_chi(this,m)
            
            m = this.sim.constrain_parameters(m);

            if isempty(this.d_current) || isempty(this.m_current) || ...
                                          (norm(m(:)-this.m_current(:)) > 1e-13)
                this.d_current = this.sim.generate_data(m);
                this.sim.display_data(this.d_current);
                this.m_current = m;
            end
            
            % Misfit
            grad_M = this.misfit.grad_M(this.d_current);
            grad_M = grad_M/numel(this.d_current);
            
            % Regularization
            grad_R = zeros(size(grad_M));
            for ri = 1:length(this.regularization)
                grad_R = grad_R + this.regularization{ri}.dR(this.m_current);
            end
            grad_R = grad_R/numel(this.m_current);
            
            this.grad_M_current = grad_M;
            this.grad_R_current = grad_R;
            
            grad_chi = grad_M + grad_R;
        end
        
        % ---------------------------- History ----------------------------
        
        function initialize_history(this)
            
            % Initialize history
            this.history.chi          = zeros(this.max_iters+1,1);
            this.history.M            = zeros(this.max_iters+1,1);
            this.history.R            = zeros(this.max_iters+1,1);
            this.history.grad_chi     = cell(this.max_iters+1,1);
            this.history.grad_M       = cell(this.max_iters+1,1);
            this.history.grad_R       = cell(this.max_iters+1,1);
            this.history.angle_chi_R  = zeros(this.max_iters+1,1);
            this.history.angle_chi_M  = zeros(this.max_iters+1,1);
            this.history.angle_M_R    = zeros(this.max_iters+1,1);
            this.history.m            = cell(this.max_iters+1,1);
            this.history.dm           = cell(this.max_iters+1,1);
            this.history.search_angle = zeros(this.max_iters+1,1);
        end
        
        function save_history(this,k,m,dm)

            this.history.chi(k+1) = this.M_current + this.R_current;
            this.history.M(k+1)   = this.M_current;
            this.history.R(k+1)   = this.R_current;
            
            grad_M   = this.grad_M_current;
            grad_R   = this.grad_R_current;
            grad_chi = grad_M + grad_R;
            
            this.history.grad_chi{k+1}  = grad_chi;
            this.history.grad_M{k+1}    = grad_M;
            this.history.grad_R{k+1}    = grad_R;

            this.history.angle_chi_R(k+1) = real(acos((grad_R'*grad_chi)/(norm(grad_R)*norm(grad_chi))));
            this.history.angle_chi_M(k+1) = real(acos((grad_M'*grad_chi)/(norm(grad_M)*norm(grad_chi))));
            this.history.angle_M_R(k+1)   = real(acos((grad_M'*grad_R)/(norm(grad_M)*norm(grad_R))));

            this.history.m{k+1}            = m;
            this.history.dm{k+1}           = dm;
            this.history.search_angle(k+1) = real(acos((-dm(:)'*this.history.grad_chi{k})/(norm(dm(:))*norm(this.history.grad_chi{k}))));
        end

        % ------------------------- Visualization -------------------------

        function setup_figures(this)

            if this.display.history.on
                this.prepare_plots('history');
            end
            if this.display.line_search.on
                this.prepare_plots('line search');
            end
        end
        
        function prepare_plots(this,option)
            
            if strcmpi(option,'history')
                this.display.history.fig = figure('name','Minimization History','NumberTitle','off');
                fh_pos = get(this.display.history.fig,'Position');
                fh_pos(1) = fh_pos(1)-0.5*fh_pos(3);
                fh_pos(3) = 2*fh_pos(3);
                set(this.display.history.fig,'Position',fh_pos);

                ax{1} = subplot(3,2,1);
                ax{2} = subplot(3,2,2);
                ax{3} = subplot(3,2,3);
                ax{4} = subplot(3,2,4);
                ax{5} = subplot(3,2,5);
                ax{6} = subplot(3,2,6);

                title(ax{1},'Objective (\chi), misfit (M) and regularization (R) function values','Fontsize',10)
                xlabel(ax{1},'k')
                xlim(ax{1},[0 1.3*(this.max_iters+1)])

                plot(ax{2},1:this.max_iters,zeros(this.max_iters,1),'k--')
                hold(ax{2},'on');
                plot(ax{2},1:this.max_iters,10*ones(this.max_iters,1),'k:')
                plot(ax{2},1:this.max_iters,-10*ones(this.max_iters,1),'k:')
                hold(ax{2},'off');
                title(ax{2},'Change of \chi, M and R relative to previous iteration','Fontsize',10)
                xlabel(ax{2},'k')
                ylabel(ax{2},'% change')
                yticks = [-20 -10 0 10 20];
                ytick_labels = {'-20','-10','0','10','20'};
                set(ax{2},'YTick',yticks,'YTickLabel',ytick_labels);
                xlim(ax{2},[0 1.3*(this.max_iters+1)])
                ylim(ax{2},[1.1*yticks(1) 1.1*yticks(end)]);

                plot(ax{3},1:this.max_iters+1,log10(this.cutoff_tol)*ones(this.max_iters+1,1),'k--')
                title(ax{3},'Norms of gradients of objective, misfit and regularization functions','Fontsize',10)
                xlabel(ax{3},'k')
                xlim(ax{3},[0 1.3*(this.max_iters+1)])
                yticks = -10:1;
                ytick_labels{1} = ' ';
                for i = 2:length(yticks)-2
                    if mod(i,2) == 0
                        ytick_labels{i} = ' ';
                    else
                        ytick_labels{i} = ['1e' num2str(i-length(yticks)+1)];
                    end
                end
                if mod(length(yticks),2) == 0
                    ytick_labels{length(yticks)-1} = '1';
                    ytick_labels{length(yticks)}   = ' ';
                else
                    ytick_labels{length(yticks)-1} = ' ';
                    ytick_labels{length(yticks)}   = '10';
                end
                set(ax{3},'YTick',yticks,'YTickLabel',ytick_labels);
                ylim(ax{3},[1.1*yticks(1) 1.1*yticks(end)]);

                plot(ax{4},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
                title(ax{4},'Angles \theta between gradients of objective, misfit and regularization functions','Fontsize',10)
                xlabel(ax{4},'k')
                ylabel(ax{4},'\theta (in rad)')
                xlim(ax{4},[0 1.3*(this.max_iters+1)])
                ylim(ax{4},[-0.15 pi+0.15])
                set(ax{4},'YTick',[0 pi/2 pi]);
                set(ax{4},'YTickLabel',{'0','pi/2','pi'});

                plot(ax{5},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
                hold(ax{5},'on')
                plot(ax{5},1:this.max_iters+1,0*ones(this.max_iters+1,1),'k--')
                hold(ax{5},'off')
                title(ax{5},'Angles \theta between steepest descent (-\nabla\chi) and search (\Delta m) directions','Fontsize',10)
                xlabel(ax{5},'k')
                ylabel(ax{5},'\theta (in rad)')
                xlim(ax{5},[0 1.3*(this.max_iters+1)])
                ylim(ax{5},[-0.15 pi/2+0.15])
                set(ax{5},'YTick',[0 pi/4 pi/2]);
                set(ax{5},'YTickLabel',{'0','pi/4','pi/2'});
                
                plot(ax{6},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
                title(ax{6},'Norms of steepest descent (-\nabla\chi) and search (\Delta m) directions','Fontsize',10)
                xlabel(ax{6},'k')
                xlim(ax{6},[0 1.3*(this.max_iters+1)])
                yticks = -4:1;
                for i = 1:length(yticks)-2
                    if mod(i,2) == 0
                        ytick_labels{i} = ' ';
                    else
                        ytick_labels{i} = ['1e' num2str(i-length(yticks)+1)];
                    end
                end
                if mod(length(yticks),2) == 0
                    ytick_labels{length(yticks)-1} = '1';
                    ytick_labels{length(yticks)}   = ' ';
                else
                    ytick_labels{length(yticks)-1} = ' ';
                    ytick_labels{length(yticks)}   = '10';
                end
                set(ax{6},'YTick',yticks,'YTickLabel',ytick_labels);
                ylim(ax{6},[1.1*yticks(1) 1.1*yticks(end)]);

                this.display.history.axis{1} = ax{1};
                this.display.history.axis{2} = ax{2};
                this.display.history.axis{3} = ax{3};
                this.display.history.axis{4} = ax{4};
                this.display.history.axis{5} = ax{5};
                this.display.history.axis{6} = ax{6};
            elseif strcmpi(option,'line search')
                this.display.line_search.fig  = figure('Name','Line Search','NumberTitle','off');
                this.display.line_search.axis = axes;
                this.line_search.set('axes handle',this.display.line_search.axis);
            end
            
            drawnow
        end
        
        function plot_history(this)
            
            for i = 1:6, ax{i} = this.display.history.axis{i}; end

            p1_R   = plot(ax{1},1:this.max_iters+1,this.history.R,'g');
            hold(ax{1},'on');
            p1_M   = plot(ax{1},1:this.max_iters+1,this.history.M,'r');
            p1_chi = plot(ax{1},1:this.max_iters+1,this.history.chi,'b');
            hold(ax{1},'off');
            title(ax{1},'Objective (\chi), misfit (M) and regularization (R) function values','Fontsize',10)
            xlabel(ax{1},'k')
            xlim(ax{1},[0 1.3*(this.max_iters+1)])
            legend(ax{1},[p1_chi p1_M p1_R],'\chi','M','R')

            chi_rel = (this.history.chi(2:end) - this.history.chi(1:end-1))./this.history.chi(1:end-1);
            chi_rel(isnan(chi_rel)) = 0;
            chi_rel(isinf(chi_rel)) = 0;

            M_rel = (this.history.M(2:end) - this.history.M(1:end-1))./this.history.M(1:end-1);
            M_rel(isnan(M_rel)) = 0;
            M_rel(isinf(M_rel)) = 0;

            R_rel = (this.history.R(2:end) - this.history.R(1:end-1))./this.history.R(1:end-1);
            R_rel(isnan(R_rel)) = 0;
            R_rel(isinf(R_rel)) = 0;

            plot(ax{2},1:this.max_iters,zeros(this.max_iters,1),'k--')
            hold(ax{2},'on');
            p2_R   = plot(ax{2},1:this.max_iters,100*R_rel,'g');
            p2_M   = plot(ax{2},1:this.max_iters,100*M_rel,'r');
            p2_chi = plot(ax{2},1:this.max_iters,100*chi_rel,'b');
            plot(ax{2},1:this.max_iters,10*ones(this.max_iters,1),'k:')
            plot(ax{2},1:this.max_iters,-10*ones(this.max_iters,1),'k:')
            hold(ax{2},'off');
            title(ax{2},'Change of \chi, M and R relative to previous iteration','Fontsize',10)
            xlabel(ax{2},'k')
            ylabel(ax{2},'% change')
            yticks = [-20 -10 0 10 20];
            ytick_labels = {'-20','-10','0','10','20'};
            set(ax{2},'YTick',yticks,'YTickLabel',ytick_labels);
            xlim(ax{2},[0 1.3*(this.max_iters+1)])
            ylim(ax{2},[1.1*yticks(1) 1.1*yticks(end)]);
            legend(ax{2},[p2_chi p2_M p2_R],'\chi','M','R')

            grad_R_norms   = zeros(this.max_iters+1,1);
            grad_M_norms   = zeros(this.max_iters+1,1);
            grad_chi_norms = zeros(this.max_iters+1,1);
            for i = 1:this.max_iters+1
                grad_R_norms(i)   = norm(this.history.grad_R{i});
                grad_M_norms(i)   = norm(this.history.grad_M{i});
                grad_chi_norms(i) = norm(this.history.grad_chi{i});
            end
            plot(ax{3},1:this.max_iters+1,log10(this.cutoff_tol)*ones(this.max_iters+1,1),'k--')
            hold(ax{3},'on');
            p3_R   = plot(ax{3},1:this.max_iters+1,log10(grad_R_norms),'g');
            p3_M   = plot(ax{3},1:this.max_iters+1,log10(grad_M_norms),'r');
            p3_chi = plot(ax{3},1:this.max_iters+1,log10(grad_chi_norms),'b');
            hold(ax{3},'off');
            title(ax{3},'Norms of gradients of objective, misfit and regularization functions','Fontsize',10)
            xlabel(ax{3},'k')
            xlim(ax{3},[0 1.3*(this.max_iters+1)])
            yticks = sort([log10(this.cutoff_tol) -10:1]); [~,ind] = find(yticks == log10(this.cutoff_tol)); yticks = yticks(ind:end);
            ytick_labels{1} = 'cutoff';
            for i = 2:length(yticks)-2
                if mod(i,2) == 0
                    ytick_labels{i} = ' ';
                else
                    ytick_labels{i} = ['1e' num2str(i-length(yticks)+1)];
                end
            end
            if mod(length(yticks),2) == 0
                ytick_labels{length(yticks)-1} = '1';
                ytick_labels{length(yticks)}   = ' ';
            else
                ytick_labels{length(yticks)-1} = ' ';
                ytick_labels{length(yticks)}   = '10';
            end
            set(ax{3},'YTick',yticks,'YTickLabel',ytick_labels);
            ylim(ax{3},[1.1*yticks(1) 1.1*yticks(end)]);
            legend(ax{3},[p3_chi p3_M p3_R],'||\nabla\chi||','||\nabla M||','||\nabla R||')

            plot(ax{4},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
            hold(ax{4},'on');
            p4_3 = plot(ax{4},1:this.max_iters+1,this.history.angle_M_R,'g');
            p4_2 = plot(ax{4},1:this.max_iters+1,this.history.angle_chi_R,'r');
            p4_1 = plot(ax{4},1:this.max_iters+1,this.history.angle_chi_M,'b');
            hold(ax{4},'off');
            title(ax{4},'Angles \theta between gradients of objective, misfit and regularization functions','Fontsize',10)
            xlabel(ax{4},'k')
            ylabel(ax{4},'\theta (in rad)')
            xlim(ax{4},[0 1.3*(this.max_iters+1)])
            ylim(ax{4},[0 pi])
            set(ax{4},'YTick',[0 pi/2 pi]);
            set(ax{4},'YTickLabel',{'0','pi/2','pi'});
            legend(ax{4},[p4_1 p4_2 p4_3],'\theta(\nabla\chi,\nabla M)','\theta(\nabla\chi,\nabla R)','\theta(\nabla M,\nabla R)')

            plot(ax{5},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
            hold(ax{5},'on');
            plot(ax{5},1:this.max_iters+1,0*ones(this.max_iters+1,1),'k--')
            p5 = plot(ax{5},1:this.max_iters+1,this.history.search_angle,'r');
            hold(ax{5},'off');
            title(ax{5},'Angles \theta between steepest descent (-\nabla\chi) and search (\Delta m) directions','Fontsize',10)
            xlabel(ax{5},'k')
            ylabel(ax{5},'\theta (in rad)')
            xlim(ax{5},[0 1.3*(this.max_iters+1)])
            ylim(ax{5},[-0.15 pi/2+0.15])
            set(ax{5},'YTick',[0 pi/4 pi/2]);
            set(ax{5},'YTickLabel',{'0','pi/4','pi/2'});
            
            dm_norms       = zeros(this.max_iters+1,1);
            grad_chi_norms = zeros(this.max_iters+1,1);
            for i = 1:this.max_iters+1
                dm_norms(i)       = norm(this.history.dm{i}(:));
                grad_chi_norms(i) = norm(this.history.grad_chi{i});
            end
            plot(ax{6},1:this.max_iters+1,pi/2*ones(this.max_iters+1,1),'k--')
            hold(ax{6},'on');
            p6_chi = plot(ax{6},1:this.max_iters+1,log10(grad_chi_norms),'b');
            p6_dm  = plot(ax{6},1:this.max_iters+1,log10(dm_norms),'r');
            hold(ax{6},'off');
            title(ax{6},'Norms of steepest descent (-\nabla\chi) and search (\Delta m) directions','Fontsize',10)
            xlabel(ax{6},'k')
            xlim(ax{6},[0 1.3*(this.max_iters+1)])
            yticks = -4:1;
            for i = 1:length(yticks)-2
                if mod(i,2) == 0
                    ytick_labels{i} = ' ';
                else
                    ytick_labels{i} = ['1e' num2str(i-length(yticks)+1)];
                end
            end
            if mod(length(yticks),2) == 0
                ytick_labels{length(yticks)-1} = '1';
                ytick_labels{length(yticks)}   = ' ';
            else
                ytick_labels{length(yticks)-1} = ' ';
                ytick_labels{length(yticks)}   = '10';
            end
            set(ax{6},'YTick',yticks,'YTickLabel',ytick_labels);
            ylim(ax{6},[1.1*yticks(1) 1.1*yticks(end)]);
            legend(ax{6},[p6_chi p6_dm],'||\nabla\chi||','||\Delta m||')
            
            drawnow
        end
 
        function switch_plot(this,option,on)
            
            if strcmpi(option,'history')
                if this.display.history.on
                    if on
                        set(0,'currentfigure',this.display.history.fig{1});
                    else
                        this.display.on = 0;
                        if isfield(this.display.history.fig)
                            for i = 1:length(this.display.history.fig)
                                close(this.display.history.fig{i});
                            end
                            this.display.history.fig = []; this.display.history.axis = [];
                        end
                    end
                else
                    if on
                        this.display.history.on = 1;
                        this.prepare_plots('history');
                    end
                end
            elseif strcmpi(option,'line search')
                if this.display.line_search.on
                    if on
                        set(0,'currentfigure',this.display.line_search.fig{1});
                    else
                        this.display.on = 0;
                        if isfield(this.display.line_search.fig)
                            for i = 1:length(this.display.line_search.fig)
                                close(this.display.line_search.fig{i});
                            end
                            this.display.line_search.fig = []; this.display.line_search.axis = [];
                            this.line_search.set('axis handle',this.line_search.axis);
                        end
                    end
                else
                    if on
                        this.display.line_search.on = 1;
                        this.prepare_plots('line search');
                    end
                end
            end
        end
    end
end
