classdef Integrator < handle

    properties (Access = public)
        name = '';          % Name of integrator
        
        t = [];             % Integration times
        max_dt = 0;         % Maximum time step
        
        y = [];             % Solution of forward problem
        
        PDE = [];           % Link to PDE to be integrated
        sim = [];           % Link to Simulation object
        
        data_times = [];    % Times at which data is collected; integration
                            % times must include these times
        
        visualization = 0;  % Solution is displayed during the
                            % integration procedure if true
    end
    
    methods (Abstract)
        
        y      = get_forward_solution(this,varargin);
        lambda = get_adjoint_solution(this);
        
        % -----------------------------------------------------------------

        yy = solve_forward(this,varargin);
        lambda = solve_adjoint(this,theta);
        
        % -----------------------------------------------------------------

        Ty = compute_T(this,y);
        Ty = compute_T_inv(this,varargin);      % Same as 'solve_forward'

        % ======================= First Derivatives =======================
        
        % ------------------- Derivative wrt solution ---------------------

        q      = compute_DtDy(this,w,varargin);
        theta  = compute_DtDy_transpose(this,lambda,varargin);
        w      = compute_DtDy_inv(this,q);
        lambda = compute_DtDy_transpose_inv(this,theta,varargin);
        
        DtDmw  = compute_DtDm(this,w,varargin);
        DtDmTw = compute_DtDm_transpose(this,w,varargin);
        
        DtDy0w  = compute_DtDy0(this,w,varargin);
        DtDy0Tw = compute_DtDy0_transpose(this,w,varargin);
        
        DyDm  = compute_DyDm(this,w);
        DyDmt = compute_DyDm_transpose(this,w);
        
        % ======================= Second Derivatives ======================
        
        q = compute_DDtDyDy(this,v,w);
        q = compute_DDtDyDy_transpose(this,v,w);
        q = compute_DDtDmDy(this,v,w);
        q = compute_DDtDmDy_transpose(this,v,w);
        q = compute_DDtDy0Dy(this,v,w);
        q = compute_DDtDy0Dy_transpose(this,v,w);
        
        q = compute_DDtDyDm(this,v,w);
        q = compute_DDtDyDm_transpose(this,v,w);
        q = compute_DDtDmDm(this,v,w);
        q = compute_DDtDmDm_transpose(this,v,w);
        q = compute_DDtDy0Dm(this,v,w);
        q = compute_DDtDy0Dm_transpose(this,v,w);
        
        q = compute_DDtDyDy0(this,v,w);
        q = compute_DDtDyDy0_transpose(this,v,w);
        q = compute_DDtDmDy0(this,v,w);
        q = compute_DDtDmDy0_transpose(this,v,w);
        q = compute_DDtDy0Dy0(this,v,w);
        q = compute_DDtDy0Dy0_transpose(this,v,w);
        
        % =================================================================
        
    end
    
    methods (Access = public)

        % ----------------------- Mutator functions -----------------------

        function set(this, varargin)
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('Options for ''set'':');
                    disp('  ''PDE''');
                    disp('  ''integration times''');
                    disp('  ''max dt'' or '' max time step''');
                    disp('  ''visualization'' or ''display''');
                elseif strcmpi(parameter,'sim')
                    this.sim = varargin{i+1};
                elseif strcmpi(parameter,'PDE')
                    this.set_PDE(varargin{i+1});
                elseif strcmpi(parameter,'integration times')
                    this.set_int_times(varargin{i+1});
                elseif strcmpi(parameter,'max dt') || strcmpi(parameter,'max time step')
                    this.set_max_dt(varargin{i+1});
                end
            end
        end
        
        function set_PDE(this,PDE)
            this.PDE = PDE;
        end

        function set_int_times(this,varargin)

            data_t = this.data_times;
            if nargin > 1
                data_t = varargin{1}(:);
            end
            if isempty(this.data_times), this.data_times = data_t; end
            
            dt = this.max_dt;
            if nargin > 2
                dt = varargin{1};
            end
            
            data_t  = [0; data_t];
            tt = [];
            for ti = 1:length(data_t)-1
                Dt = data_t(ti+1) - data_t(ti);
                N  = ceil((Dt-sqrt(eps))/dt);
                tt = [tt; linspace(data_t(ti),data_t(ti+1),N+1)'];
                tt = tt(1:end-1);
            end
            this.t  = [tt; data_t(end)];
            
            this.sim.compute_Q();
        end
        
        function set_max_dt(this,max_dt)
            this.max_dt = max_dt;
            this.set_int_times();
        end

        % ---------------------- Accessor functions -----------------------

        function out = get(this,parameter,varargin)
            out = [];
            
            if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                disp('Options for ''get'':');
                disp('  ''PDE''');
                disp('  ''integration times''');
                disp('  ''num time steps'' or ''number of time steps'' or ''K''');
                disp('  ''max dt'' or ''max time step''');
                disp('  ''solution'' or ''forward solution''');
                disp('  ''adjoint solution''');
                disp('  ''visualization'' or ''display''');
            elseif strcmpi(parameter,'PDE')
                out = this.get_PDE();
            elseif strcmpi(parameter,'integration times')
                out = this.get_time_steps();
            elseif strcmpi(parameter,'num time steps') || strcmpi(parameter,'K') ...
                                        || strcmpi(parameter, 'number of time steps')
                out = length(this.t)-1;
            elseif strcmpi(parameter,'max dt') || strcmpi(parameter,'max time step')
                out = this.max_dt();
            elseif strcmpi(parameter,'solution') || strcmpi(parameter,'forward solution')
                out = this.get_forward_solution(varargin{:});
            elseif strcmpi(parameter,'adjoint solution')
                out = this.get_adjoint_solution(varargin{:});
            elseif strcmpi(parameter,'visualization') || strcmpi(parameter,'display')
                out = this.visualization();
            end
        end
        
        function PDE = get_PDE(this)
            PDE = this.PDE;
        end
        
        function t = get_time_steps(this)
            t = this.t;
        end
        
        function ti = get_time_indices(this,t)
            [~,ti] = ismember(t,this.t);
            ti = ti-1;
        end

        % -----------------------------------------------------------------
        
        function Dt = compute_Dt(this,w,varargin)
            
            yy = [];
            supplied_solution = 0;
            transpose = 0;
            inverse = 0;
            derivative = '';
            for i = 1:2:nargin-2
                parameter = varargin{i};
                if strcmpi(parameter,'y')
                    supplied_solution = 1;
                    yy = varargin{i+1};
                elseif strcmpi(parameter, 'transpose')
                    transpose = varargin{i+1};
                elseif strcmpi(parameter, 'inverse')
                    inverse = varargin{i+1};
                elseif strcmpi(parameter, 'derivative')
                    derivative = varargin{i+1};
                end
            end
            
            if strcmpi(derivative,'dy')
                if supplied_solution
                    if transpose
                        if inverse
                            Dt = this.compute_DtDy_transpose_sol_inv(w,yy);
                        else
                            Dt = this.compute_DtDy_transpose_sol(w,yy);
                        end
                    else
                        if inverse
                            Dt = this.compute_DtDy_sol_inv(w,yy);
                        else
                            Dt = this.compute_DtDy_sol(w,yy);
                        end
                    end
                else
                    if transpose
                        if inverse
                            Dt = this.compute_DtDy_transpose_inv(w);
                        else
                            Dt = this.compute_DtDy_transpose(w);
                        end
                    else
                        if inverse
                            Dt = this.compute_DtDy_inv(w);
                        else
                            Dt = this.compute_DtDy(w);
                        end
                    end
                end
            elseif strcmpi(derivative,'dm')
                if supplied_solution
                    disp('HELP')
                    if transpose
                        Dt = this.compute_DtDm_transpose_sol(w,yy);
                    else
                        Dt = this.compute_DtDm_sol(w,yy);
                    end
                else
                    if transpose
                        Dt = this.compute_DtDm_transpose(w);
                    else
                        Dt = this.compute_DtDm(w);
                    end
                end
            end
        end

        function DyDmt = compute_DyDmt(this,w)
            DyDmt = this.compute_DyDm_transpose(w);
        end

        % -------------------------- Diagnostics --------------------------
        
        function results = run_tests(this,varargin)
            
            if nargin > 1
                s = varargin{1};
            else
                s.derivative.DtDy    = 0;   s.inverse.T      = 0;
                s.derivative.DtDm    = 0;   s.inverse.DtDy   = 0;
                s.derivative.DyDm    = 0;   s.inverse.DtDy_t = 0;

                s.transpose.DtDy     = 0;
                s.transpose.DtDy_inv = 0;   % This is the adjoint test
                s.transpose.DtDm     = 0; 
                s.transpose.DyDm     = 0;
            end
            
            pass.derivative.DtDy    = 0;   pass.inverse.T      = 0;
            pass.derivative.DtDm    = 0;   pass.inverse.DtDy   = 0;
            pass.derivative.DyDm    = 0;   pass.inverse.DtDy_t = 0;
            
            pass.transpose.DtDy     = 0;
            pass.transpose.DtDy_inv = 0;
            pass.transpose.DtDm     = 0;
            pass.transpose.DyDm     = 0;
            
            N  = this.PDE.get('N');
            Nm = this.PDE.get('Nm');
            K  = length(this.t)-1;
            
            test_derivative = s.derivative.DtDy + s.derivative.DtDy_inv + s.derivative.DtDm + s.derivative.DyDm;
            test_transpose  = s.transpose.DtDy + s.transpose.DtDy_inv + ...
                              s.transpose.DtDm + s.transpose.DyDm;
            test_inverse    = s.inverse.T + s.inverse.DtDy + s.inverse.DtDy_t;
            
            q = 0.05*rand(N,K,this.s+1);
            y = this.compute_T_inv(q);
                        
%             Y1 = reshape(y,N,K,this.s+1);
%             for k = 1:K
%                 for ii = 1:5
%                     Y{k}{ii} = Y1(:,k,ii);
%                 end
%             end
%             Q1 = reshape(q,N,K,this.s+1);
%             for k = 1:K
%                 for ii = 1:5
%                     Q{k}{ii} = Q1(:,k,ii);
%                 end
%             end
            
            if test_derivative
                fprintf('Testing derivatives of T:');
                
                Ty    = q;
                dy    = 0.1*norm(this.y(:))*rand(N,K,this.s+1);
                dq    = 0.005*rand(N,K,this.s+1);
                dm    = this.PDE.get('dm');
                old_m = this.PDE.get('m');

                if s.derivative.DtDy,     DtDy_times_dY     = this.compute_DtDy(dy); end
                if s.derivative.DtDy_inv, DtDy_inv_times_dQ = this.compute_DtDy_inv(dq); end
                if s.derivative.DtDm,     DtDm_times_dM     = this.compute_DtDm(dm); end
                if s.derivative.DyDm,     DyDm_times_dM     = this.compute_DyDm(dm); end

                for i = 1:10
                    i
                    hy = 1*0.1^(i+1);
                    hq = 1*0.1^(i+1);
                    hm = 1*0.1^(i+1);
                    
                    if s.derivative.DtDy
                        q_pert = this.compute_T(y+hy*dy);
                        
                        err_Dt(i)   = norm(q_pert(:)-q(:))/norm(q(:));
                        err_DtDy(i) = norm(q_pert(:)-q(:)-hy*DtDy_times_dY(:))/norm(q(:));
                    end
                    
                    if s.derivative.DtDy_inv
                        y_pert = this.compute_T_inv(q+hq*dq);
                        this.compute_T_inv(q);

                        err_Dt_inv(i)   = norm(y_pert(:)-y(:))/norm(y(:));
                        err_DtDy_inv(i) = norm(y_pert(:)-y(:)-hq*DtDy_inv_times_dQ(:))/norm(y(:));
                    end
                    
                    if s.derivative.DtDm
                        this.PDE.set('m',old_m + hm*dm);
                        q_pert = this.compute_T(y);
                        this.PDE.set('m',old_m);

                        err_Dt(i)   = norm(q_pert(:)-q(:))/norm(q(:));
                        err_DtDm(i) = norm(q_pert(:)-q(:)-hm*DtDm_times_dM(:))/norm(q(:));
                    end
                    
                    if s.derivative.DyDm
                        this.PDE.set('m',old_m + hm*dm);
                        y_pert = this.compute_T_inv(q);
                        this.PDE.set('m',old_m);
                    
                        err_dY(i)   = norm(y_pert(:)- y(:))/norm(y(:));
                        err_DyDm(i) = norm(y_pert(:)- y(:)-hm*DyDm_times_dM(:))/norm(y(:));
                    end
                end
                
                if s.derivative.DtDy
                    err_Dt = err_Dt(:); err_Dt
                    err_DtDy = err_DtDy(:); err_DtDy
                end
                
                if s.derivative.DtDy_inv
                    err_Dt_inv = err_Dt_inv(:); err_Dt_inv
                    err_DtDy_inv = err_DtDy_inv(:); err_DtDy_inv
                end
                
                if s.derivative.DtDm
                    err_Dt = err_Dt(:); err_Dt
                    err_DtDm = err_DtDm(:); err_DtDm
                end
                
                if s.derivative.DyDm
                    err_dY = err_dY(:); err_dY
                    err_DyDm = err_DyDm(:); err_DyDm
                end
            end
            
            if test_transpose
                fprintf('Testing transposes of t:');
                
                if s.transpose.DtDy
                    w = rand(N,K,this.s+1);
                    v = rand(N,K,this.s+1);

                    if this.PDE.Fourier_basis
                        Nx = this.PDE.grid.x.N;
                        Ny = this.PDE.grid.y.N;
                        for k = 1:K
                            for si = 1:this.s+1
                                w_temp = fft2(reshape(w(:,k,si),Nx,Ny));
                                w(:,k,si) = w_temp(:);
                                
                                v_temp = fft2(reshape(v(:,k,si),Nx,Ny));
                                v(:,k,si) = v_temp(:);
                            end
                        end
                    end

                    vt = this.compute_Dt(w,'derivative','dy');
                    wt = this.compute_Dt(v,'derivative','dy','transpose',1);
                    
                    err_DtDy_t = (v(:)'*vt(:) - w(:)'*wt(:))/abs(v(:)'*vt(:))
                end
                
                if s.transpose.DtDy_inv
                    w = rand(N,K,this.s+1);
                    v = rand(N,K,this.s+1);

                    if this.PDE.Fourier_basis
                        Nx = this.PDE.grid.x.N;
                        Ny = this.PDE.grid.y.N;
                        for k = 1:K
                            for si = 1:this.s+1
                                w_temp = fft2(reshape(w(:,k,si),Nx,Ny));
                                w(:,k,si) = w_temp(:);
                                
                                v_temp = fft2(reshape(v(:,k,si),Nx,Ny));
                                v(:,k,si) = v_temp(:);
                            end
                        end
                    end
                    
                    vt = this.compute_Dt(w,'derivative','dy','inverse',1);
                    wt = this.compute_Dt(v,'derivative','dy','inverse',1,'transpose',1);
                    
                    err_DtDy_inv_t = (vt(:)'*v(:) - wt(:)'*w(:))/abs(v(:)'*vt(:))
                end
                
                if s.transpose.DtDm
                    w = rand(Nm,1);
                    v = rand(N,K,this.s+1);
                    
                    if this.PDE.Fourier_basis
                        Nx = this.PDE.grid.x.N;
                        Ny = this.PDE.grid.y.N;
                        for k = 1:K
                            for si = 1:this.s+1
                                v_temp = fft2(reshape(v(:,k,si),Nx,Ny));
                                v(:,k,si) = v_temp(:);
                            end
                        end
                    end

                    vt = this.compute_Dt(w,'derivative','dm');
                    wt = this.compute_Dt(v,'derivative','dm','transpose',1);

                    err_DtDm_t = (vt(:)'*v(:) - wt(:)'*w(:))/abs(v(:)'*vt(:))
                end
                
                if s.transpose.DyDm
                    w = rand(Nm,1);
                    v = rand(N,K,this.s+1);
                    
                    if this.PDE.Fourier_basis
                        Nx = this.PDE.grid.x.N;
                        Ny = this.PDE.grid.y.N;
                        for k = 1:K
                            for si = 1:this.s+1
                                v_temp = fft2(reshape(v(:,k,si),Nx,Ny));
                                v(:,k,si) = v_temp(:);
                            end
                        end
                    end

                    vt = this.compute_DyDm(w);
                    wt = this.compute_DyDm_transpose(v);

                    err_DyDm_t = (vt(:)'*v(:) - wt(:)'*w(:))/abs(v(:)'*vt(:))
                end
            end
            
            if test_inverse
                fprintf('Testing inverses of T:');
                
                if s.inverse.T
                    Ty = this.compute_T(y);
                    
                    err_T_inv = abs(norm(q(:) - Ty(:))/norm(q(:)))
                end
                
                if s.inverse.DtDy
                    w = rand(N,K,this.s+1);
                    
                    theta  = this.compute_Dt(w,    'derivative','dy');
                    lambda = this.compute_Dt(theta,'derivative','dy','inverse',1);

                    err_DtDy_inv = abs(norm(w(:) - lambda(:))/norm(w(:)))
                end
                
                if s.inverse.DtDy_t
                    w = rand(N,K,this.s+1);

                    theta  = this.compute_Dt(w,    'derivative','dy','transpose',1);
                    lambda = this.compute_Dt(theta,'derivative','dy','transpose',1,'inverse',1);

                    err_DtDyt_inv = abs(norm(w(:) - lambda(:))/norm(w(:)))
                end
            end
            
            results = pass;
        end
      
    end
end
