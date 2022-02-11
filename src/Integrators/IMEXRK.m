classdef IMEXRK < Integrator

    properties (Access = public)

        scheme = 0;         % RK scheme (classical 4th-order will be used by default)
        c = [];             % Nodes of RK scheme
        s = 0;              % Number of stages of RK scheme
        
        lambda;             % Solution of adjoint problem
        
        explicit            = 1;
        diagonally_implicit = 0;
        implicit            = 0;
        
        i1 = 0
        i2 = 0
        
        Y = [];
        
        % Temporary storage
        u05 = []; u1 = [];
        R05 = []; R1 = [];
       
    end
    
    methods (Access = public)
        
        function this = IMEXRK(varargin)
            
            this.name = 'IMEXRK';
            
            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
            else
                s = settings('integrator');
            end

            if nargin > 1
                this.set(varargin{:});
            end
            
            if isempty(this.contour)
                this.setup_contour(s.contour.M,s.contour.shape);
            end
            
            if this.scheme == 0
                this.set_scheme(s.scheme);
            end
            
            if this.max_dt <= 0
                this.max_dt = s.max_dt;
            end
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this,varargin)
            
            set@Integrator(this,varargin{:});
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('  ''scheme''');
                elseif strcmpi(parameter,'scheme')
                    this.set_scheme(varargin{i+1});
                end
            end

            this.precompute_coefficients();
        end

        function set_scheme(this,scheme)
            if strcmpi(scheme, 'Euler')
                this.scheme = 1;
                this.c = 0;
            elseif strcmpi(scheme, 'Cox-Matthews')
                this.scheme = 2;
                this.c = [0; 0.5; 0.5; 1];
            elseif strcmpi(scheme, 'Krogstad')
                this.scheme = 3;
                this.c = [0; 0.5; 0.5; 1];
            elseif strcmpi(scheme, 'Hochbruck-Ostermann')
                this.scheme = 4;
                this.c = [0; 0.5; 0.5; 1; 0.5];
            end
            
            this.s = length(this.c);
            this.precompute_coefficients();
        end
               
        % ---------------------- Accessor functions -----------------------
        
        function out = get(this,parameter,varargin)

            out = get@Integrator(this,parameter,varargin{:});

            if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                disp('Options for ''get'':');
                disp('  ''scheme''');
                disp('  ''s'' or ''num stages''');
            elseif strcmpi(parameter,'scheme')
                out = this.scheme;
            elseif strcmpi(parameter,'s') || strcmpi(parameter,'num stages')
                out = this.s;
            end
        end
        
        function yt = get_solution(this,varargin)
            yt = this.get_forward_solution(varargin{:});
        end

        function yt = get_forward_solution(this,varargin)
            if nargin > 1
                yt = varargin{1}(:,:,this.s+1);
            else
                yt = this.y(:,:,this.s);
            end
            if this.PDE.L_diagonal()
                yt = this.PDE.change_basis(yt,'inverse');
            end
        end
        
        function lambda = get_adjoint_solution(this,varargin)
            if nargin > 1
                lambda = varargin{1}(:,:,this.s+1);
            else
                lambda = this.lambda(:,:,this.s);
            end
            if this.PDE.L_diagonal()
                lambda = this.PDE.change_basis(lambda,'inverse');
            end
        end

        function scheme = get_scheme(this)
            scheme = this.scheme;
        end
        
        % -----------------------------------------------------------------

        function yy = solve_forward(this,varargin)
            if nargout > 1
                yy = this.compute_T_inv(varargin{:});
            else
                this.compute_T_inv(varargin{:});
            end
        end

        function lambda = solve_adjoint(this,theta)
            lambda = this.compute_DtDy_transpose_inv(theta);
        end
        
        % -----------------------------------------------------------------

        function q = compute_T(this,y)
            
            N  = this.PDE.get('num unknowns');
            K  = length(this.t)-1;
            y0 = this.PDE.get('initial condition');

            q      = zeros(N,K,this.s+1);
            y_temp = zeros(N,this.s);
            
            this.prepare_resolvents(this.t(2) - this.t(1));

            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                if k == 1, yk = y0(:); else yk = y(:,k-1,this.s+1); end
                
                this.solve_resolvents(yk,dtk,yk);
                
                y_temp(:,1) = yk;
                for si = 2:this.s
                    y_temp(:,si) = this.exp_times(this.c(si));
                end
                q(:,k,this.s+1) = y(:,k,this.s+1) - this.exp_times(1);
                
                for si = 1:this.s
                    this.solve_resolvents(y(:,k,si),dtk,yk);
                    for sj = si+1:this.s
                        y_temp(:,sj) = y_temp(:,sj) + dtk*this.a_times(sj,si);
                    end
                    q(:,k,si)       = y(:,k,si) - this.PDE.N(y_temp(:,si),this.t(k)+this.c(si)*dtk,yk);
                    q(:,k,this.s+1) = q(:,k,this.s+1) - dtk*this.b_times(si);
                end
            end
        end
        
        function Ty = compute_T_inv(this,varargin)
            
            ah = [];
            q  = [];
            if nargin == 3
                q = varargin{1};
                ah = varargin{2};
            elseif nargin == 2
                if ishandle(varargin{1})
                    ah = varargin{1};
                else
                    q = varargin{1};
                end
            end
            
            N = this.PDE.get('N');
            K = length(this.t)-1;
            
            this.y = zeros(N,K,this.s);
            this.Y = zeros(N,K,this.s);
            
            if nargout > 0
                Ty = zeros(N,K,this.s+1);
            end
            
            y0 = this.PDE.get('initial condition');
            yk = y0(:);
            
            yc = zeros(N,this.s);
            Y  = zeros(N,this.s);

            for k = 1:K
                dtk = this.t(k+1) - this.t(k);

                for si = 1:this.s
                    Y(:,si) = this.PDE.F(yc(:,si),this.t(k)+this.c(si)*dtk,yk);
                    if ~isempty(q)
                        Y(:,si) = Y(:,si) + q(:,k,si);
                    end
                    
                    yk1 = yk1 + dtk*this.b_times(si);
                    
                    for sj = si+1:this.s
                        yc(:,sj) = yc(:,sj) + dtk*this.a_times(sj,si);
                    end
                end
                if ~isempty(q), yk1 = yk1 + q(:,k,this.s+1); end
                yk = yk1;
                
                % Store solution
                for si = 2:this.s, this.y(:,k,si-1) = yc(:,si); end
                for si = 1:this.s, this.Y(:,k,si) = Y(:,si); end
                this.y(:,k,this.s) = yk;
                
                if nargout > 0
                    for si = 1:this.s, Ty(:,k,si) = Y(:,si); end
                    Ty(:,k,this.s+1) = yk;
                end

                if ~isempty(ah)
                    this.PDE.display_solution(yk,this.t(k+1),ah);
                end
            end
        end
        
        % ======================= First Derivatives =======================
        
        % ------------------- Derivative wrt solution ---------------------

        function q = compute_DtDy(this,w,varargin)

            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            q = zeros(N,K,this.s+1);
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end

            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                
                for si = 1:s
                    if k == 1
                        yk = y0(:);
                        wk = 0*yk;
                    else
                        yk = y(:,k-1,this.s+1);
                        wk = w(:,k-1,this.s+1);
                    end
                    
                    tkc = this.t(k)+this.c(si)*dtk;

                    for sj = 1:si-1
                        wk = wk + dtk*this.a(si,sj)*w(:,k,si);
                        yk = yk + dtk*this.a(si,sj)*y(:,k,si);
                    end
                    
                    q(:,k,si) = q(:,k,si) - this.PDE.DFDy(wk,yk,tkc,yk);
                end
                
                if k == 1
                    q(:,k,this.s+1) = w(:,k,this.s+1);
                else
                    q(:,k,this.s+1) = w(:,k,this.s+1) - w(:,k-1,this.s+1);
                end
                
                for si = 1:this.s
                    q(:,k,this.s+1) = q(:,k,this.s+1) - dtk*this.b(si)*w(:,k,si);
                end
            end
        end
        
        function w = compute_DtDy_inv(this,q)

            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            w_temp = zeros(N,this.s);
            w = zeros(N,K,this.s+1);
            
            this.prepare_resolvents(this.t(2) - this.t(1),1);
            
            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                if k == 1
                    yk = y0(:);
                else
                    yk = this.y(:,k-1,this.s);
                end
                
                w_temp = 0*w_temp;
                for si = 1:this.s
                    w(:,k,si) = q(:,k,si);
                end
                w(:,k,this.s+1) = q(:,k,this.s+1);
                
                if k > 1
                    this.solve_resolvents(w(:,k-1,this.s+1),dtk,yk);
                    
                    w_temp(:,1) = w(:,k-1,this.s+1);
                    for si = 2:this.s
                        w_temp(:,si) = this.exp_times(this.c(si));
                    end
                    w(:,k,this.s+1) = w(:,k,this.s+1) + this.exp_times(1);
                end
                    
                for si = 1:this.s
                    if si == 1
                        w(:,k,1) = w(:,k,1) + this.PDE.DNDy(w_temp(:,1),yk,this.t(k),yk);
                        if this.PDE.rosenbrock && k > 1
                            w(:,k,1) = w(:,k,1) + this.PDE.DNDyk(w(:,k-1,this.s+1),yk,this.t(k),yk);
                        end
                    else
                        tkc = this.t(k) + this.c(si)*dtk;
                        w(:,k,si) = w(:,k,si) + this.PDE.DNDy(w_temp(:,si),this.y(:,k,si-1),tkc,yk);
                        if this.PDE.rosenbrock && k > 1
                            wk = w(:,k-1,this.s+1);
                            temp = this.DexpDyk_times(this.c(si),wk,yk,dtk,yk);
                            for sj = 1:si-1
                                temp = temp + dtk*this.DaDyk_times(si,sj,wk,this.Y(:,k,sj),dtk,yk);
                            end
                            w(:,k,si) = w(:,k,si) + this.PDE.DNDy(temp,this.y(:,k,si-1),tkc,yk);
                            w(:,k,si) = w(:,k,si) + this.PDE.DNDyk(w(:,k-1,this.s+1),this.y(:,k,si-1),tkc,yk);
                        end
                    end
                    this.solve_resolvents(w(:,k,si),dtk,yk);
                    for sj = si+1:this.s
                        w_temp(:,sj) = w_temp(:,sj) + dtk*this.a_times(sj,si);
                    end
                    w(:,k,this.s+1) = w(:,k,this.s+1) + dtk*this.b_times(si);
                end
                
                if this.PDE.rosenbrock && k > 1
                    wk = w(:,k-1,this.s+1);
                    w(:,k,this.s+1) = w(:,k,this.s+1) + this.DexpDyk_times(1,wk,yk,dtk,yk);
                    for si = 1:this.s
                        w(:,k,this.s+1) = w(:,k,this.s+1) + dtk*this.DbDyk_times(si,wk,this.Y(:,k,si),dtk,yk);
                    end
                end
            end
        end
        
        function theta = compute_DtDy_transpose(this,lambda,varargin)

            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end
            
            theta = zeros(N,K,this.s+1);
            theta(:,K,this.s+1) = lambda(:,K,this.s+1);
               
            for k = K:-1:1
                dtk = this.t(k+1) - this.t(k);
                
                if k < K
                    theta(:,k,this.s+1) = theta(:,k,this.s+1) - lambda(:,k+1,this.s+1);
                    for si = 1:s
                        if k == 1
                            yk = y0(:);
                        else
                            yk = y(:,k-1,this.s+1);
                        end

                        tkc = this.t(k)+this.c(si)*dtk;

                        for sj = 1:si-1
                            yk = yk + dtk*this.a(si,sj)*y(:,k,sj);
                        end
                        theta(:,k,this.s+1) = theta(:,k,this.s+1) - this.PDE.DFDy(lambda(:,k+1,si),yk,tkc,yk,1);
                    end
                end
                
                for si = 1:s
                    theta(:,k,si) = lambda(:,k,si) - tkc*this.b(si)*lambda(:,k,this.s+1);
                    for sj = 1:s
                        if k == 1
                            yk = y0(:);
                        else
                            yk = y(:,k-1,this.s+1);
                        end

                        tkc = this.t(k)+this.c(sj)*dtk;

                        for sk = 1:sj-1
                            yk = yk + dtk*this.a(sj,sk)*y(:,k,sk);
                        end
                        theta(:,k,si) = theta(:,k,si) - dtk*this.a(sj,si)*this.PDE.DFDy(lambda(:,k,si),yk,tkc,yk,1);
                    end
                end
            end
        end
        
        function lambda = compute_DtDy_transpose_inv(this,theta,varargin)

            ah = [];
            if nargin > 2
                ah = varargin{1};
            end
            
            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            full_theta = 0;
            if size(theta,3) > 1
                full_theta = 1;
            else
                theta = reshape(theta,[],K);
            end
            
            lambda = zeros(N,K,this.s+1);

            for k = K:-1:1
                dtk = this.t(k+1) - this.t(k);
                
                if k == 1
                    yk = y0(:);
                else
                    yk = this.y(:,k-1,this.s);
                end

                if full_theta
                    lambda(:,k,this.s+1) = lambda(:,k,this.s+1) + theta(:,k,this.s+1);
                else
                    lambda(:,k,this.s+1) = lambda(:,k,this.s+1) + theta(:,k);
                end
                this.solve_resolvents(lambda(:,k,this.s+1),dtk,yk,1);
                
                for si = this.s:-1:1
                    lambda(:,k,si) = dtk*this.b_times(si);
                    if full_theta
                        lambda(:,k,si) = lambda(:,k,si) + theta(:,k,si);
                    end
                end
                
                if k > 1
                    l_temp = this.exp_times(1);
                    if this.PDE.rosenbrock
                        l_temp = l_temp + this.DexpDyk_times(1,lambda(:,k,this.s+1),yk,dtk,yk,1);
                        for si = 1:this.s
                            l_temp = l_temp + dtk*this.DbDyk_times(si,lambda(:,k,this.s+1),this.Y(:,k,si),dtk,yk,1);
                        end
                    end
                    lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + l_temp;
                end

                for si = this.s:-1:2
                    tkc = this.t(k) + this.c(si)*dtk;
                    l_temp = this.PDE.DNDy(lambda(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                    this.solve_resolvents(l_temp,dtk,yk,1);

                    for sj = si-1:-1:1
                        lambda(:,k,sj) = lambda(:,k,sj) + dtk*this.a_times(si,sj);
                    end

                    if k > 1
                        lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + this.exp_times(this.c(si));
                        if this.PDE.rosenbrock
                            lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + this.PDE.DNDyk(lambda(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                            lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + this.DexpDyk_times(this.c(si),l_temp,yk,dtk,yk,1);
                            for sj = si-1:-1:1
                                lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + dtk*this.DaDyk_times(si,sj,l_temp,this.Y(:,k,sj),dtk,yk,1);
                            end
                        end
                    end
                end
                
                if k > 1
                    lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + this.PDE.DNDy(lambda(:,k,1),yk,this.t(k),yk,1);
                    
                    if this.PDE.rosenbrock
                        lambda(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) + this.PDE.DNDyk(lambda(:,k,1),yk,this.t(k),yk,1);
                    end
                end

                if ~isempty(ah)
                    this.PDE.display_solution(lambda(:,k,this.s+1),this.t(k+1),ah);
                end
            end
        end
        
        % ---------------- Derivative wrt model parameters ----------------
        
        function DtDmw = compute_DtDm(this,w,varargin)
            
            N  = this.PDE.get('N');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            DtDmw = zeros(N,K,this.s+1);
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end

            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                for si = 1:this.s
                    if k == 1
                        yk = y0(:);
                    else
                        yk = y(:,k-1,this.s+1);
                    end

                    tkc = this.t(k)+this.c(si)*dtk;

                    for sj = 1:si-1
                        yk = yk + dtk*this.a(si,sj)*y(:,k,sj);
                    end
                    DtDmw(:,k,si) = -this.PDE.DFDm(w,yk,tkc,yk);
                end
            end
        end

        function DtDmTw = compute_DtDm_transpose(this,w,varargin)
            
            Nm = this.PDE.get('Nm');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            DtDmTw = zeros(Nm,1);
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end
                
            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                for si = 1:this.s
                    if k == 1
                        yk = y0(:);
                    else
                        yk = y(:,k-1,this.s+1);
                    end

                    tkc = this.t(k)+this.c(si)*dtk;

                    for sj = 1:si-1
                        yk = yk + dtk*this.a(si,sj)*y(:,k,sj);
                    end
                    DtDmTw = DtDmTw - this.PDE.DFDm(w,yk,tkc,yk,1);
                end
            end
        end
        
        % ---------------- Derivative wrt initial condition ---------------
        
        function DtDy0w = compute_DtDy0(this,w,varargin)
            
            N  = this.PDE.get('N');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            DtDy0w = zeros(N,K,this.s+1);
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end
            
            for si = 1:this.s
                yk = y0(:);

                tkc = this.t(1)+this.c(si)*(this.t(2) - this.t(1));

                for sj = 1:si-1
                    yk = yk + dtk*this.a(si,sj)*y(:,1,sj);
                end
                
                DtDy0w(:,1,si) = -this.PDE.DFDy0(w,yk,tkc,yk);
            end
            
            DtDy0w(:,1,this.s+1) = -w;
        end

        function DtDy0Tw = compute_DtDy0_transpose(this,w,varargin)
            
            N  = this.PDE.get('N');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            DtDy0Tw = -w;
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end
                
            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                for si = 1:this.s
                    if k == 1
                        yk = y0(:);
                    else
                        yk = y(:,k-1,this.s+1);
                    end

                    tkc = this.t(k)+this.c(si)*dtk;

                    for sj = 1:si-1
                        yk = yk + dtk*this.a(si,sj)*y(:,k,sj);
                    end
                    DtDy0Tw = DtDy0Tw - this.PDE.DFDy0(w,yk,tkc,yk,1);
                end
            end
        end
        
        % -----------------------------------------------------------------

        function DyDm = compute_DyDm(this,w)
            DyDm = -1*this.compute_DtDy_inv(this.compute_DtDm(w));
        end
        
        function DyDmt = compute_DyDm_transpose(this,w)
            DyDmt  = -1*this.compute_DtDm_transpose(this.solve_adjoint(w));
        end
        
        % ======================= Second Derivatives ======================
        
        function q = compute_DDtDyDy(this,k,v,w)
            
            N  = this.PDE.get('N');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            q = zeros(N,this.s+1);
            
            if nargin > 2
                y = varargin{1};
            else
                y = this.y;
            end
            
            if k == 1
                yk = y0;
                wk = 0*yk;
                vk = 0*yk;
            else
                yk = y0;
                wk = 0*yk;
                vk = 0*yk;
            end
            
            for si = 1:this.s
                q(:,si) = ;
            end
        end
        
        function q = compute_DDtDyDy_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDmDy(this,v,w)
            
        end
        
        function q = compute_DDtDmDy_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dy(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dy_transpose(this,v,w)
            
        end
        
        % -----------------------------------------------------------------

        function q = compute_DDtDyDm(this,v,w)
            
        end
        
        function q = compute_DDtDyDm_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDmDm(this,v,w)
            
        end
        
        function q = compute_DDtDmDm_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dm(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dm_transpose(this,v,w)
            
        end
        
        % -----------------------------------------------------------------

        function q = compute_DDtDyDy0(this,v,w)
            
        end
        
        function q = compute_DDtDyDy0_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDmDy0(this,v,w)
            
        end
        
        function q = compute_DDtDmDy0_transpose(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dy0(this,v,w)
            
        end
        
        function q = compute_DDtDy0Dy0_transpose(this,v,w)
            
        end
    end
end
