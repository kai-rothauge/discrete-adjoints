classdef Simulation < handle
    
    properties (SetAccess = public)
        
        integrator   = [];          % 'Integrator' object
        PDE          = [];          % 'PDE' object
        minimization = [];          % 'Minimization' object
        
        observations = [];          % Observations (data, times and noise level)
        Q = [];
            
        % Display settings
        display = [];
    end
    
    methods (Access = public)
        function this = Simulation(varargin)
            
            if nargin == 1 && isnumeric(varargin{1})
                this.PDE = SwiftHohenberg(varargin{1});
            end

            this.integrator   = ETDRK();
            this.minimization = Minimization();

            this.observations = settings('observations');
            if nargin > 1
                this.set(varargin{:});
            end
                
            if isempty(this.PDE)
                this.PDE          = SwiftHohenberg();
            end
            if isempty(this.integrator)
                this.integrator   = ETDRK1();
            end
            if isempty(this.minimization)
                this.minimization = Minimization();
            end
            if isempty(this.observations)
                this.observations = settings('observations');
            end

            this.integrator.set('sim',this,'PDE',this.PDE,'integration times',this.observations.times);
            this.minimization.set('sim',this);
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this,varargin)
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('Options for ''set'':');
                    disp('  ''integrator'' or ''time integrator''');
                    disp('  ''PDE''');
                    disp('  ''minimization''');
                    disp('  ''observation times'' or ''obs times''');
                    disp('  ''parameters'' or ''m''');
                elseif strcmpi(parameter,'integrator') || strcmpi(parameter,'time integrator')
                    this.set_integrator(varargin{i+1});
                elseif strcmpi(parameter,'PDE')
                    this.set_PDE(varargin{i+1});
                elseif strcmpi(parameter,'minimization')
                    this.set_minimization(varargin{i+1});
                elseif strcmpi(parameter,'observation times') || strcmpi(parameter,'obs times')
                    this.set_obs_times(varargin{i+1});
                elseif strcmpi(parameter,'parameters') || strcmpi(parameter,'m')
                    this.PDE.set('m',varargin{i+1});
                end
            end
        end
        
        function set_obs_times(this,obs_times)
            this.observations.times = obs_times;
            this.integrator.set('integration times',obs_times);
        end
        
        % ---------------------- Accessor functions -----------------------
        
        function out = get(this,parameter)
            
            out = this.PDE.get(parameter);
 
            if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                disp('  ''num data points'' or ''Nd'' or ''number of data points''');
            elseif strcmpi(parameter,'num data points') || strcmpi(parameter,'Nd') ...
                                        || strcmpi(parameter,'number of data points')
                out = length(this.observations.data(:));
            end
        end
        
        % -----------------------------------------------------------------

        function G = get_parameter_gradient(this)
            G = this.PDE.get_parameter_gradient();
        end
        
        function compute_Q(this)
            ti = this.integrator.get_time_indices(this.observations.times);
            K  = this.integrator.get('K');
            N  = this.PDE.get('N');

            I_N = speye(N);
            I_K = speye(K); I_K = I_K(ti,:);
            
            this.Q = kron(I_K,I_N);
        end
        
        % ---------------------------- Display ----------------------------

        function switch_plot(this,option,on)
            
            if strcmpi(option,'IC')
                if this.display.IC.on
                    if on, set(0,'currentfigure',this.display.IC.fig);
                    else
                        close(this.display.IC.fig);
                        this.display.IC.on = 0;
                        this.display.IC.fig  = [];
                        this.display.IC.axis = [];
                    end
                else
                    if on, this.display.IC.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'forward solution') || strcmpi(option,'solution') || strcmpi(option,'y')
                if this.display.forward_solution.on
                    if on, set(0,'currentfigure',this.display.forward_solution.fig);
                    else
                        close(this.display.forward_solution.fig);
                        this.display.forward_solution.on = 0;
                        this.display.forward_solution.fig  = [];
                        this.display.forward_solution.axis = [];
                    end
                else
                    if on, this.display.forward_solution.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'adjoint solution') || strcmpi(option,'adjoint') || strcmpi(option,'lambda')
                if this.display.adjoint_solution.on
                    if on, set(0,'currentfigure',this.display.adjoint_solution.fig);
                    else
                        close(this.display.adjoint_solution.fig);
                        this.display.adjoint_solution.on = 0;
                        this.display.adjoint_solution.fig  = [];
                        this.display.adjoint_solution.axis = [];
                    end
                else
                    if on, this.display.adjoint_solution.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'actual parameters') || strcmpi(option,'parameters actual')
                if this.display.parameters.actual.on
                    if on, set(0,'currentfigure',this.display.parameters.actual.fig{1});
                    else
                        this.display.parameters.actual.on = 0;
                        for i = 1:length(this.display.parameters.actual.fig)
                            close(this.display.parameters.actual.fig{i});
                            this.display.parameters.actual.fig{i}  = [];
                            this.display.parameters.actual.axis{i} = [];
                        end
                    end
                else
                    if on, this.display.parameters.actual.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'initial parameters') || strcmpi(option,'parameters initial')
                if this.display.parameters.initial.on
                    if on, set(0,'currentfigure',this.display.parameters.initial.fig{1});
                    else
                        this.display.parameters.initial.on = 0;
                        for i = 1:length(this.display.parameters.initial.fig)
                            close(this.display.parameters.initial.fig{i});
                            this.display.parameters.initial.fig{i}  = [];
                            this.display.parameters.initial.axis{i} = [];
                        end
                    end
                else
                    if on, this.display.parameters.initial.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'current parameters') || strcmpi(option,'parameters current')
                if this.display.parameters.current.on
                    if on, set(0,'currentfigure',this.display.parameters.current.fig{1});
                    else
                        this.display.parameters.current.on = 0;
                        for i = 1:length(this.display.parameters.current.fig)
                            close(this.display.parameters.current.fig{i});
                            this.display.parameters.current.fig{i}  = [];
                            this.display.parameters.current.axis{i} = [];
                        end
                    end
                else
                    if on, this.display.parameters.current.on = 1; this.prepare_plot(option); end
                end
            elseif strcmpi(option,'data actual')     || strcmpi(option,'actual data')  || ...
                   strcmpi(option,'data initial')    || strcmpi(option,'initial data') || ...
                   strcmpi(option,'data current')    || strcmpi(option,'current data') || ...
                   strcmpi(option,'data difference') || strcmpi(option,'difference data')
                if this.display.data.actual.on  || this.display.data.initial.on || ...
                   this.display.data.current.on || this.display.data.difference.on
                    if on, set(0,'currentfigure',this.display.data.fig);
                    else
                        close(this.display.data.fig);
                        this.display.data.fig             = [];
                        this.display.data.actual.on       = 0;
                        this.display.data.initial.on      = 0;
                        this.display.data.current.on      = 0;
                        this.display.data.difference.on   = 0;
                        this.display.data.actual.axis     = [];
                        this.display.data.initial.axis    = [];
                        this.display.data.current.axis    = [];
                        this.display.data.difference.axis = [];
                    end
                else
                    if on
                        if strcmpi(option,'data actual') || strcmpi(option,'actual data')
                            this.display.data.actual.on = 1;
                        elseif strcmpi(option,'data initial') || strcmpi(option,'initial data')
                            this.display.data.initial.on = 1;
                        elseif strcmpi(option,'data current') || strcmpi(option,'current data')
                            this.display.data.current.on = 1;
                        elseif strcmpi(option,'data difference') || strcmpi(option,'difference data')
                            this.display.data.difference.on = 1;
                        end
                        this.prepare_plot(option);
                    end
                end
            end
        end
        
        function prepare_plot(this,option)
            
            if strcmpi(option,'IC')
                this.display.IC.fig  = figure('Name','Initial Condition','NumberTitle','off');
                this.display.IC.axis = axes;
            elseif strcmpi(option,'forward solution') || strcmpi(option,'solution') || strcmpi(option,'y')
                this.display.forward_solution.fig  = figure('Name','Forward Solution','NumberTitle','off');
                this.display.forward_solution.axis = axes;
            elseif strcmpi(option,'adjoint solution') || strcmpi(option,'adjoint') || strcmpi(option,'lambda')
                this.display.adjoint_solution.fig  = figure('Name','Adjoint Solution','NumberTitle','off');
                this.display.adjoint_solution.axis = axes;
            elseif strcmpi(option,'actual parameters') || strcmpi(option,'parameters actual')
                m_names = this.PDE.get_parameter_names();
                for j = 1:length(m_names)
                    fig_name = ['Parameter ' m_names{j}.name];
                    this.display.parameters.actual.fig{j}  = figure('Name',fig_name,'NumberTitle','off');
                    this.display.parameters.actual.axis{j} = axes;
                end
            elseif strcmpi(option,'initial parameters') || strcmpi(option,'parameters initial')
                m_names = this.PDE.get_parameter_names();
                
                for j = 1:length(m_names)
                    fig_name = ['Parameter ' m_names{j}.name ' (Initial Guess)'];
                    this.display.parameters.initial.fig{j}  = figure('Name',fig_name,'NumberTitle','off');
                    this.display.parameters.initial.axis{j} = axes;
                end
            elseif strcmpi(option,'current parameters') || strcmpi(option,'parameters current')
                m_names = this.PDE.get_parameter_names();
                
                for j = 1:length(m_names)
                    fig_name = ['Parameter ' m_names{j}.name ' (Current Estimate)'];
                    this.display.parameters.current.fig{j}  = figure('Name',fig_name,'NumberTitle','off');
                    this.display.parameters.current.axis{j} = axes;
                end
            elseif strcmpi(option,'data')
                this.display.data.fig = figure('Name','Data','NumberTitle','off');
                if this.display.data.actual.on && ~this.display.data.initial.on && ...
                   this.display.data.current.on && ~this.display.data.difference.on
                    this.display.data.axis.actual = axes;
                else
                    this.display.data.axis.actual     = subplot(2,2,1);
                    this.display.data.axis.initial    = subplot(2,2,2);
                    this.display.data.axis.current    = subplot(2,2,3);
                    this.display.data.axis.difference = subplot(2,2,4);
                end
            end
        end
        
        function prepare_plots(this,s)

            this.display.IC.on   = 0;
            this.display.IC.fig  = [];
            this.display.IC.axis = [];
            
            this.display.forward_solution.on   = 0;
            this.display.forward_solution.fig  = [];
            this.display.forward_solution.axis = [];
            
            this.display.adjoint_solution.on   = 0;
            this.display.adjoint_solution.fig  = [];
            this.display.adjoint_solution.axis = [];
                
            m_names = this.PDE.get_parameter_names();

            this.display.parameters.actual.on  = 0;
            this.display.parameters.initial.on = 0;
            this.display.parameters.current.on = 0;

            this.display.parameters.actual.fig  = [];
            this.display.parameters.initial.fig = [];
            this.display.parameters.current.fig = [];
            for j = 1:length(m_names)
                this.display.parameters.actual.fig{j}   = [];
                this.display.parameters.initial.fig{j}  = [];
                this.display.parameters.current.fig{j}  = [];
            end

            this.display.parameters.actual.axis   = [];
            this.display.parameters.initial.axis  = [];
            this.display.parameters.current.axis  = [];
            for j = 1:length(m_names)
                this.display.parameters.actual.axis{j}  = [];
                this.display.parameters.initial.axis{j} = [];
                this.display.parameters.current.axis{j} = [];
            end
            
            this.display.data.fig             = []; 
            
            this.display.data.actual.on       = 0; 
            this.display.data.actual.axis     = []; 
            
            this.display.data.initial.on      = 0; 
            this.display.data.initial.axis    = []; 
            
            this.display.data.current.on      = 0; 
            this.display.data.current.axis    = []; 
            
            this.display.data.difference.on   = 0; 
            this.display.data.difference.axis = []; 

            this.display.minimization.history.on = 0;
            this.display.minimization.line_search.on = 0;
  
            if isfield(s,'forward')
                
                on = 0;
                if isfield(s.forward,'IC'), on = s.forward.IC; end
                this.switch_plot('IC',on);
                
                on = 0;
                if isfield(s.forward,'solution'), on = s.forward.solution; end
                this.switch_plot('forward solution',on);
                
                on = 0;
                if isfield(s.forward,'parameters'), on = s.forward.parameters; end
                this.switch_plot('parameters actual',on);
                
            elseif isfield(s,'inverse')
                
                if isfield(s.inverse,'parameters')
                    on = 0;
                    if isfield(s.inverse.parameters,'actual'), on = s.inverse.parameters.actual; end
                    this.switch_plot('parameters actual',on);

                    on = 0;
                    if isfield(s.inverse.parameters,'initial'), on = s.inverse.parameters.initial; end
                    this.switch_plot('parameters initial',on);

                    on = 0;
                    if isfield(s.inverse.parameters,'current'), on = s.inverse.parameters.current; end
                    this.switch_plot('parameters current',on);
                end
                    
                if isfield(s.inverse,'solution')
                    on = 0;
                    if isfield(s.inverse.solution,'forward'), on = s.inverse.solution.forward; end
                    this.switch_plot('forward solution',on);
                    
                    on = 0;
                    if isfield(s.inverse.solution,'adjoint'), on = s.inverse.solution.adjoint; end
                    this.switch_plot('adjoint solution',on);
                end
                
                if isfield(s.inverse,'minimization')
                    on = 0;
                    if isfield(s.inverse.minimization,'history'), on = s.inverse.minimization.history; end
                    this.minimization.set('display history',on);
                    
                    on = 0;
                    if isfield(s.inverse.minimization,'line_search'), on = s.inverse.minimization.line_search; end
                    this.minimization.set('display line search',on);
                end

                on = 0;
                if isfield(s.inverse,'data'), on = s.inverse.data; end
                this.switch_plot('data',on);
            end
        end
        
        function results = run(this,option)
            if option == 1
                this.prepare_plots(settings('display','forward'));
                
                fprintf('Solving forward problem ............................... ');
                t_start = cputime;
                this.PDE.set('parameters','actual');
                this.PDE.display_parameters(this.display.parameters.actual.axis);
                this.PDE.display_initial_condition(this.display.IC.axis);
                if ~isempty(this.display.forward_solution.fig)
                    figure(this.display.forward_solution.fig);
                end
                this.integrator.solve_forward(this.display.forward_solution.axis);
                results = this.integrator.get('solution');
                fprintf('done (%3.4fs)\n',cputime - t_start);
            elseif option == 2
                this.prepare_plots(settings('display','inverse'));
                
                fprintf('-----------------------------------------------------------------------\n');
                fprintf('                     Recovering parameters\n');
                fprintf('-----------------------------------------------------------------------\n');
                fprintf('Generating data ....................................... ');
                t_start = cputime;
                this.generate_observations();
                fprintf('done (%3.4fs)\n',cputime - t_start);

                this.PDE.set('parameters','initial guess');
                this.PDE.display_parameters(this.display.parameters.initial.axis);
                this.minimization.set('m_ref',this.PDE.get('parameters'),...
                                      'd_obs',this.observations.data);

                results = this.minimization.run();
            elseif option == 3
                fprintf('-----------------------------------------------------------------------\n');
                fprintf('                     Comparing orders of accuracy\n');
                fprintf('-----------------------------------------------------------------------\n');

                max_dt = this.integrator.get('max time step');
                
                dts = max_dt./[8 4 2 1];
                schemes = {'Euler','Cox-Matthews','Krogstad','Hochbruck-Ostermann'};
                N = 1;
                
                E_y      = zeros(length(schemes),length(dts)-1,N);
                E_lambda = zeros(length(schemes),length(dts)-1,N);
                E_grad   = zeros(length(schemes),length(dts)-1,N);

                header_1 = '    -----------';
                header_2 = '              |';
                header_3 = '              |';
                header_4 = '    -----------';
                footer   = '    -----------';
                for si = 1:length(schemes)
                    if strcmpi(schemes{si},'Euler')
                        header_1 = [header_1 '-------------'];
                        header_2 = [header_2 '             '];
                        header_3 = [header_3 '    Euler    '];
                        header_4 = [header_4 '-------------'];
                        footer   = [footer   '-------------'];
                    elseif strcmpi(schemes{si},'Cox-Matthews')
                        header_1 = [header_1 '-------------'];
                        header_2 = [header_2 '     Cox-    '];
                        header_3 = [header_3 '   Matthews  '];
                        header_4 = [header_4 '-------------'];
                        footer   = [footer   '-------------'];
                    elseif strcmpi(schemes{si},'Krogstad')
                        header_1 = [header_1 '-------------'];
                        header_2 = [header_2 '             '];
                        header_3 = [header_3 '   Krogstad  '];
                        header_4 = [header_4 '-------------'];
                        footer   = [footer   '-------------'];
                    elseif strcmpi(schemes{si},'Hochbruck-Ostermann')
                        header_1 = [header_1 '-------------'];
                        header_2 = [header_2 '  Hochbruck- '];
                        header_3 = [header_3 '  Ostermann  '];
                        header_4 = [header_4 '-------------'];
                        footer   = [footer   '-------------'];
                    end
                end
                header = [header_1 '-\n' header_2 '|\n' header_3 '|\n' header_4 '-\n'];
                footer = [footer '-\n\n'];
                
                for n = 1:N
                    this.PDE.set_IC('random');
                    this.integrator.set('scheme','Krogstad');
                    this.integrator.set('max dt',0.5*dts(1));
                
                    fprintf('\n                              Run %d of %d\n\n',n,N);
                    fprintf('Generating data ....................................... ');
                    t_start = cputime;
                    this.PDE.set('parameters','actual');
                    this.integrator.solve_forward();
                    d_obs = this.extract_data();
                    fprintf('done (%3.4fs)\n',cputime - t_start);

                    fprintf('Solving "exact" forward problem ....................... ');
                    t_start = cputime;
                    this.PDE.set('parameters','initial guess');
                    this.integrator.solve_forward();
                    y_exact = this.integrator.get('solution');
                    fprintf('done (%3.4fs)\n',cputime - t_start);

                    fprintf('Solving "exact" adjoint problem ....................... ');
                    t_start = cputime;
                    d_sim = this.extract_data();
                    dDdYt = this.Q'*(d_sim(:) - d_obs(:));
                    if this.PDE.L_diagonal()
                        dDdYt = reshape(dDdYt,this.PDE.get('N'),[]);
                        dDdYt = this.PDE.change_basis(dDdYt);
                        dDdYt = dDdYt(:)/this.PDE.get('N');
                    end
                    dDdYt = reshape(dDdYt,this.PDE.get('N'),[]);
                    l = this.integrator.solve_adjoint(dDdYt);
                    lambda_exact = this.integrator.get('adjoint solution',l);
                    fprintf('done (%3.4fs)\n',cputime - t_start);

                    fprintf('Computing "exact" gradient ............................ ');
                    t_start = cputime;
                    grad_exact = this.integrator.compute_DyDmt(dDdYt);
                    fprintf('done (%3.4fs)\n',cputime - t_start);

                    for si = 1:length(schemes)
                        fprintf(['\nEvaluating ' schemes{si} ' scheme:\n']);
                        this.integrator.set('scheme',schemes{si});

                        for i = 1:length(dts)
                            fprintf('  Time step dt = %2.3e:\n',dts(i));
                            this.integrator.set('max dt',dts(i));

                            fprintf('    Solving forward problem ........................... ');
                            t_start = cputime;
                            this.PDE.set('parameters','initial guess');
                            this.integrator.solve_forward();
                            y{i} = this.integrator.get('solution');
                            fprintf('done (%3.4fs)\n',cputime - t_start);

                            fprintf('    Solving adjoint problem ........................... ');
                            t_start = cputime;
                            l = this.integrator.solve_adjoint(dDdYt(:,2^(i):2^(i):end));
                            lambda{i} = this.integrator.get('solution',l);
                            fprintf('done (%3.4fs)\n',cputime - t_start);

                            fprintf('    Computing gradient ................................ ');
                            t_start = cputime;
                            grad{i} = this.integrator.compute_DyDmt(dDdYt(:,2^(i):2^(i):end));
                            fprintf('done (%3.4fs)\n',cputime - t_start);
                        end

                        K = size(y{end},2);
                        I = length(dts);
                        for i = 1:I
                            err_y{i}      = y{i}(:,2^(I-i):2^(I-i):end)      - y_exact(:,2^I:2^I:end);
                            err_lambda{i} = lambda{i}(:,2^(I-i):2^(I-i):end) - lambda_exact(:,2^I:2^I:end);
                            err_grad{i}   = grad{i} - grad_exact;
                            if i > 1
                                E_y(si,i-1,n)      = log2(norm(err_y{i}(:))/norm(err_y{i-1}(:)));
                                E_lambda(si,i-1,n) = log2(norm(err_lambda{i}(:))/norm(err_lambda{i-1}(:)));
                                E_grad(si,i-1,n)   = log2(norm(err_grad{i}(:))/norm(err_grad{i-1}(:)));
                                
%                                 if (si == 1 && (E_y(si,i-1,n) < 0.5))
%                                     E_y(si,i-1,n) = -1;
%                                 elseif (si > 1 && (E_y(si,i-1,n) < 3.5))
%                                     E_y(si,i-1,n) = -1;
%                                 end
%                                 
%                                 if (si == 1 && (E_lambda(si,i-1,n) < 0.5))
%                                     E_lambda(si,i-1,n) = -1;
%                                 elseif (si > 1 && (E_lambda(si,i-1,n) < 3.5))
%                                     E_lambda(si,i-1,n) = -1;
%                                 end
%                                 
%                                 if (si == 1 && (E_grad(si,i-1,n) < 0.5))
%                                     E_grad(si,i-1,n) = -1;
%                                 elseif (si > 1 && (E_grad(si,i-1,n) < 3.5))
%                                     E_grad(si,i-1,n) = -1;
%                                 end
                            end
                        end
                        
                        fprintf('\n');
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the forward solution for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_y(i*dt)|/|err_y(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_y(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_y(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the adjoint solution for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_lambda(i*dt)|/|err_lambda(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_lambda(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_lambda(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the gradient for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_gradient(i*dt)|/|err_gradient(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_grad(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_grad(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                    end
                end
                
                for si = 1:length(schemes)
                    for i = 1:length(dts)-1
                        E_y(si,i,1)      = mean(E_y(si,i,(find(E_y(si,i,:) ~= -1))));
                        E_lambda(si,i,1) = mean(E_lambda(si,i,(find(E_lambda(si,i,:) ~= -1))));
                        E_grad(si,i,1)   = mean(E_grad(si,i,(find(E_grad(si,i,:) ~= -1))));
                    end
                end

                fprintf('\n');
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the forward solution for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_y(i*dt)|/|err_y(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_y(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_y(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the adjoint solution for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_lambda(i*dt)|/|err_lambda(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_lambda(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_lambda(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                fprintf('Approximate orders of accuracy of the gradient for selected schemes:\n');
                fprintf('        E(i,j) = log2(|err_gradient(i*dt)|/|err_gradient(j*dt)|\n');
                fprintf(header);
                for i = 1:length(dts)-1
                    fprintf('     E(%2d,%2d) |',2^(i),2^(i-1));
                    for si = 1:length(schemes)
                        if isnan(E_grad(si,i,1))
                            fprintf('      *      ');
                        else
                            fprintf('    %3.4f   ',E_grad(si,i,1));
                        end
                    end
                    fprintf('|\n');
                end
                fprintf(footer);
                
                keyboard
            elseif option == 5
                fprintf('Computing Jacobian J:\n'); 
                fprintf('Generating data ....................................... ');
                t_start = cputime;
                this.generate_observations();
                this.integrator.set('visualization','on');
                fprintf('done (%3.4fs)\n',cputime - t_start);
                
                K  = length(this.observations.times);
                N  = this.get('N');
                Nm = this.get('Nm');
                
                max_rows = K;
                max_cols = 1024;

                E = speye(Nm);
                
                Ir = 1:ceil(K/max_rows):max(K,max_rows);
                Ic = 1:ceil(Nm/max_cols):max(Nm,max_cols);

                J = zeros(length(Ir),length(Ic));

                for i = 1:length(Ic)
                    i
                    J_col = this.J_times(E(:,Ic(i)));
                    for k = 1:K
                        J(k,i) = norm(J_col((k-1)*N+1:k*N));
                    end

                    figure(11)
                    imagesc(J)
                    colorbar
                    xlabel('Parameters')
                    ylabel('Data')
                    drawnow
                end

                results = J;
            elseif option == 4
                disp('Running tests:'); 

                s = settings('diagnostics');

%                 results.Jacobian   = this.run_tests(this);
%                 results.PDE        = this.PDE.run_tests();
                results.integrator = this.integrator.run_tests(s.integration);
%                 results = this.minimization.misfit.run_tests();
%                 results = this.minimization.regularization{1}.run_tests();
%                 results = this.minimization.regularization{2}.run_tests();
            end
        end
        
        function Jv = J_times(this,v,varargin)
            
            transpose = 0;
            if nargin > 2
                transpose = varargin{1};
            end
            
            if transpose
                dDdYt = this.Q'*v(:);
                if this.PDE.L_diagonal()
                    dDdYt = reshape(dDdYt,this.PDE.get('N'),[]);
                    dDdYt = this.PDE.change_basis(dDdYt);
                    dDdYt = dDdYt(:)/this.PDE.get('N');
                end
                Jv = this.integrator.compute_dYdMt(dDdYt);
            else
                y  = this.integrator.get_solution(this.integrator.compute_dYdM(v));
                Jv = this.Q*y(:);
            end
        end
        
        function D = extract_data(this,varargin)
            if nargin > 1
                y = varargin{1};
            else
                y = this.integrator.get('solution');
            end
            D = reshape(this.Q*y(:),[],length(this.observations.times));
        end
        
        function D = generate_data(this,varargin)
            this.PDE.set('parameters',varargin{:});
            this.integrator.solve_forward();
            
            % Get solution at observation times
            D = this.extract_data();
        end
        
        function data = D(this,m)
            this.PDE.set('m',m);
            this.integrator.solve_forward();
            
            % Get solution at observation times
            data = this.extract_data();
        end
        
        function m = constrain_parameters(this,varargin)
            m = this.PDE.constrain_parameters(varargin{:});
        end
        
        function update_plots(this,varargin)
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'m')||strcmpi(parameter,'parameters')||strcmpi(parameter,'model parameters')
                    this.PDE.display_parameters(varargin{i+1},this.display.parameters.current.axis);
                elseif strcmpi(parameter,'data')
                    this.display_data(varargin{i+1});
                end
            end
        end
        
        function display_data(this,varargin)
            
            if ~isempty(this.display.data.current.axis)
                grid = this.PDE.get('grid');

                if nargin > 1
                    D = varargin{1};
                    if nargin > 2, di = varargin{2}; else di = size(D,2); end
                else
                    D = this.observations.data;
                    di = size(D,2);
                end

                x = linspace(grid.x.min,grid.x.max,grid.x.N+1);
                y = linspace(grid.y.min,grid.y.max,grid.y.N+1);

                imagesc(x,y,reshape(D(:,di),grid.x.N,grid.y.N),'Parent',this.display.data.current.axis);
                colormap(this.axes_handle,'parula')
                colorbar(this.axes_handle)
                xlabel(this.axes_handle,'x','Fontsize',10)
                ylabel(this.axes_handle,'y','rotation',0,'Fontsize',10)
                axis(this.axes_handle,'equal','tight')
                drawnow
            end
        end

        function generate_observations(this)
            this.PDE.set('parameters','actual');
            this.PDE.display_parameters(this.display.parameters.actual.axis);
            this.integrator.solve_forward();
            
            % Get solution at observation times
            D = this.extract_data();
            N = sqrt(this.PDE.grid.x.N*this.PDE.grid.y.N);

            % Add noise
            for i = 1:size(D,2)
                Di = D(:,i);
                D(:,i) = D(:,i) + 0.01*this.observations.noise*norm(Di)*randn(size(Di))/N;
            end
            this.observations.data = D;
            this.minimization.set('d_obs',D);
        end
        
        function result = run_tests(this,N)

            result = [];
            
            disp('Testing Jacobian:')
            
            m  = this.PDE.get('m');
            dm = this.PDE.get('dm');

            D             = this.D(m); D = D(:);
            dDdM_times_dM = this.J_times(dm);

            err_D = zeros(10,1); 
            err_J = zeros(10,1); 
            
            for i = 1:10
                h = 0.1^(i+2);
                
                D_pert = this.D(m + h*dm);
                D_pert = D_pert(:);

                D_res_1 = D_pert - D;
                D_res_2 = D_pert - D - h*dDdM_times_dM;
                err_D(i) = norm(D_res_1(:))/norm(D(:));
                err_J(i) = norm(D_res_2(:))/norm(D(:));
            end
            
            w = dm;
            v = 0.1*rand(size(D));

            vt = this.J_times(w);
            wt = this.J_times(v,1);

            err_Jt = (vt(:)'*v(:) - wt(:)'*w(:))/abs(v(:)'*vt(:))
            
            result.derivative.D = err_D;
            result.derivative.J = err_J;
            result.transpose.J  = err_Jt;
        end
    end
end
