classdef ETDRK < Integrator

    properties (Access = public)

        scheme = 0;         % ETD scheme (Cox-Matthews will be used by default)
        c = [];             % Nodes of ETD scheme
        s = 0;              % Number of stages of ETD scheme
        
        lambda;             % Solution of adjoint problem
        
        contour = [];       % Used for contour integration of phi-functions

        i1 = 0
        i2 = 0
        
        Y = [];
        
        % Temporary storage
        u05 = []; u1 = [];
        R05 = []; R1 = [];
        
        biw  = [];
        aijw = [];
        expw = [];
        
        a_coeffs = [];
        b_coeffs = [];
        exp_coeffs = [];
    end
    
    methods (Access = public)
        
        function this = ETDRK(varargin)
            
            this.name = 'ETDRK';
            
            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
            else
                s = settings('integration');
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
            
            this.prepare_resolvents(this.t(2) - this.t(1));
            
            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                
                this.solve_resolvents(yk,dtk,yk);
                    
                yc(:,1) = yk;
                for si = 2:this.s
                    yc(:,si) = this.exp_times(this.c(si));
                end
                yk1 = this.exp_times(1);
                
                for si = 1:this.s
                    Y(:,si) = this.PDE.N(yc(:,si),this.t(k)+this.c(si)*dtk,yk);
                    if ~isempty(q), Y(:,si) = Y(:,si) + q(:,k,si); end
                    this.solve_resolvents(Y(:,si),dtk,yk);
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
        
        % ------------------- Derivative wrt solution ---------------------

        function q = compute_DtDy(this,w)

            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            q = zeros(N,K,this.s+1);
            this.prepare_resolvents(this.t(2) - this.t(1));

            for k = 1:K
                dtk = this.t(k+1) - this.t(k);
                
                for si = 1:this.s
                    q(:,k,si) = w(:,k,si);
                end
                q(:,k,this.s+1) = w(:,k,this.s+1);
                
                if k == 1
                    yk = y0(:);
                else
                    wk = w(:,k-1,this.s+1);
                    yk = this.y(:,k-1,this.s);

                    this.solve_resolvents(wk,dtk,yk);
                    q(:,k,1) = q(:,k,1) - this.PDE.DNDy(wk,yk,this.t(k),yk);
                    for si = 2:this.s
                        tkc = this.t(k)+this.c(si)*dtk;
                        temp = this.exp_times(this.c(si));
                        q(:,k,si) = q(:,k,si) - this.PDE.DNDy(temp,this.y(:,k,si-1),tkc,yk);
                    end
                    q(:,k,this.s+1) = q(:,k,this.s+1) - this.exp_times(1);
                
                    if this.PDE.rosenbrock
                        q(:,k,this.s+1) = q(:,k,this.s+1) - this.DexpDyk_times(1,wk,yk,dtk,yk);
                        for si = 1:this.s
                            q(:,k,this.s+1) = q(:,k,this.s+1) - dtk*this.DbDyk_times(si,wk,this.Y(:,k,si),dtk,yk);
                        end
                        
                        q(:,k,1) = q(:,k,1) - this.PDE.DNDyk(wk,yk,this.t(k),yk);
                        for si = 2:this.s
                            tkc = this.t(k)+this.c(si)*dtk;
                            temp = this.DexpDyk_times(this.c(si),wk,yk,dtk,yk);
                            for sj = 1:si-1
                                temp = temp + dtk*this.DaDyk_times(si,sj,wk,this.Y(:,k,sj),dtk,yk);
                            end
                            q(:,k,si) = q(:,k,si) - this.PDE.DNDy(temp,this.y(:,k,si-1),tkc,yk);
                            q(:,k,si) = q(:,k,si) - this.PDE.DNDyk(wk,this.y(:,k,si-1),tkc,yk);
                        end
                    end
                end
                
                for si = 1:this.s
                    this.solve_resolvents(w(:,k,si),dtk,yk);
                    for sj = si+1:this.s
                        tkc = this.t(k)+this.c(sj)*dtk;
                        temp = this.a_times(sj,si);
                        q(:,k,sj) = q(:,k,sj) - dtk*this.PDE.DNDy(temp,this.y(:,k,sj-1),tkc,yk);
                    end
                    q(:,k,this.s+1) = q(:,k,this.s+1) - dtk*this.b_times(si);
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
        
        function theta = compute_DtDy_transpose(this,lambda)

            N  = this.PDE.get('num unknowns');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;
            
            theta = zeros(N,K,this.s+1);
            theta(:,K,this.s+1) = lambda(:,K,this.s+1);
               
            for k = K:-1:1
                dtk = this.t(k+1) - this.t(k);
                if k == 1
                    yk = y0(:);
                else
                    yk = this.y(:,k-1,this.s);
                end
                
                this.solve_resolvents(lambda(:,k,this.s+1),dtk,yk,1);
                
                for si = this.s:-1:1
                    theta(:,k,si) = lambda(:,k,si) - dtk*this.b_times(si);
                end
                if k > 1
                    theta(:,k-1,this.s+1) = lambda(:,k-1,this.s+1) - this.exp_times(1);
                    if this.PDE.rosenbrock
                        theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - this.DexpDyk_times(1,lambda(:,k,this.s+1),yk,dtk,yk,1);
                        for si = 1:this.s
                            theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - dtk*this.DbDyk_times(si,lambda(:,k,this.s+1),this.Y(:,k,si),dtk,yk,1);
                        end
                    end
                end
                
                for si = this.s:-1:1
                    tkc = this.t(k) + this.c(si)*dtk;
                    if si > 1
                        l_temp = this.PDE.DNDy(lambda(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                    else
                        l_temp = this.PDE.DNDy(lambda(:,k,si),yk,tkc,yk,1);
                    end
                    this.solve_resolvents(l_temp,dtk,yk,1);
                    
                    for sj = si-1:-1:1
                        theta(:,k,sj) = theta(:,k,sj) - dtk*this.a_times(si,sj);
                    end
                    if k > 1
                        if si > 1
                            theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - this.exp_times(this.c(si));
                            if this.PDE.rosenbrock
                                theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - this.PDE.DNDyk(lambda(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                                theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - this.DexpDyk_times(this.c(si),l_temp,yk,dtk,yk,1);
                                for sj = si-1:-1:1
                                    theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - dtk*this.DaDyk_times(si,sj,l_temp,this.Y(:,k,sj),dtk,yk,1);
                                end
                            end
                        else
                            theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - l_temp;
                            if this.PDE.rosenbrock
                                theta(:,k-1,this.s+1) = theta(:,k-1,this.s+1) - this.PDE.DNDyk(lambda(:,k,1),yk,this.t(k),yk,1);
                            end
                        end
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
                
                for k = 1:K
                    dtk = this.t(k+1) - this.t(k);
                    if k == 1
                        yk = y0(:);
                    else
                        yk = y(:,k-1,this.s+1);
                    end
                    for si = 1:this.s
                        y_t = this.exp_times(yk,dtk,this.c(si),'yk',yk);
                        for i = 1:si-1
                            y_t = y_t + dtk*this.a_times(y(:,k,i),si,i,dtk,'yk',yk);
                        end
                        tkc = this.t(k) + this.c(si)*dtk;
                        DtDmw(:,k,si) = -this.PDE.DNDm(w,y_t,tkc,yk);
                        
                        if this.PDE.rosenbrock && si > 1
                            temp = this.DexpDm_times(this.c(si),w,yk,dtk,yk);
                            for sj = 1:si-1
                                temp = temp + dtk*this.DaDm_times(si,sj,w,this.Y(:,k,sj),dtk,yk);
                            end
                            
                            DtDmw(:,k,si) = DtDmw(:,k,si) - this.PDE.DNDy(temp,this.y(:,k,si-1),tkc,yk);
                        end
                    end
                    if this.PDE.rosenbrock
                        DtDmw(:,k,this.s+1) = -this.DexpDm_times(1,w,yk,dtk,yk);
                        for si = 1:this.s
                            DtDmw(:,k,this.s+1) = DtDmw(:,k,this.s+1) - dtk*this.DbDm_times(si,w,this.Y(:,k,si),dtk,yk);
                        end
                    end
                end
            else
                for k = 1:K
                    dtk = this.t(k+1) - this.t(k);
                    if k == 1
                        yk = y0(:);
                    else
                        yk = this.y(:,k-1,this.s);
                    end
                    DtDmw(:,k,1) = -this.PDE.DNDm(w,yk,this.t(k),yk);
                    for si = 2:this.s
                        tkc = this.t(k) + this.c(si)*dtk;
                        DtDmw(:,k,si) = -this.PDE.DNDm(w,this.y(:,k,si-1),tkc,yk);
                    end
                    if this.PDE.rosenbrock
                        for si = 2:this.s
                            tkc = this.t(k) + this.c(si)*dtk;
                            temp = this.DexpDm_times(this.c(si),w,yk,dtk,yk);
                            for sj = 1:si-1
                                temp = temp + dtk*this.DaDm_times(si,sj,w,this.Y(:,k,sj),dtk,yk);
                            end
                            
                            DtDmw(:,k,si) = DtDmw(:,k,si) - this.PDE.DNDy(temp,this.y(:,k,si-1),tkc,yk);
                        end
                    
                        DtDmw(:,k,this.s+1) = -this.DexpDm_times(1,w,yk,dtk,yk);
                        for si = 1:this.s
                            DtDmw(:,k,this.s+1) = DtDmw(:,k,this.s+1) - dtk*this.DbDm_times(si,w,this.Y(:,k,si),dtk,yk);
                        end
                    end
                end
            end
        end

        function DtDmtw = compute_DtDm_transpose(this,w,varargin)
            
            Nm = this.PDE.get('Nm');
            y0 = this.PDE.get('initial condition');
            K  = length(this.t)-1;

            DtDmtw = zeros(Nm,1);
            
            if nargin > 2
                y = varargin{1};
                
                for k = 1:K
                    dtk = this.t(k+1) - this.t(k);
                    if k == 1
                        yk = y0(:);
                    else
                        yk = y(:,k-1,this.s+1);
                    end
                    if this.PDE.rosenbrock
                        DtDmtw = DtDmtw - this.DexpDm_times(1,w(:,k,this.s+1),yk,dtk,yk,1);
                    end
                    for si = 1:this.s
                        y_t = this.exp_times(yk,dtk,this.c(si),'yk',yk);
                        for i = 1:si-1
                            y_t = y_t + dtk*this.a_times(y(:,k,i),si,i,dtk,'yk',yk);
                        end
                        tkc = this.t(k) + this.c(si)*dtk;
                        if this.PDE.rosenbrock
                            DtDmtw = DtDmtw - dtk*this.DbDm_times(si,w(:,k,this.s+1),this.Y(:,k,si),dtk,yk,1);
                            
                            if si > 1
                                v = this.PDE.DNDy(w(:,k,this.s+1),yt,tkc,yk,1);
                                DtDmtw = DtDmtw - this.DexpDm_times(this.c(si),v,yk,dtk,yk,1);
                                for sj = 1:si-1
                                    DtDmtw = DtDmtw - dtk*this.DaDm_times(si,sj,v,this.y(:,k,sj),dtk,yk,1);
                                end
                            end
                        end
                        DtDmtw = DtDmtw - this.PDE.DNDm(w(:,k,si),y_t,tkc,yk,1);
                    end
                end
            else
                for k = 1:K
                    dtk = this.t(k+1) - this.t(k);
                    if k == 1
                        yk = y0(:);
                    else
                        yk = this.y(:,k-1,this.s);
                    end
                    DtDmtw = DtDmtw - this.PDE.DNDm(w(:,k,1),yk,this.t(k),yk,1);
                    for si = 2:this.s
                        tkc = this.t(k) + this.c(si)*dtk;
                        DtDmtw = DtDmtw - this.PDE.DNDm(w(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                    end
                    
                    if this.PDE.rosenbrock
                        DtDmtw = DtDmtw - this.DexpDm_times(1,w(:,k,this.s+1),yk,dtk,yk,1);
                        for si = 1:this.s
                            tkc = this.t(k) + this.c(si)*dtk;
                            DtDmtw = DtDmtw - dtk*this.DbDm_times(si,w(:,k,this.s+1),this.Y(:,k,si),dtk,yk,1);
                            
                            if si > 1
                                v = this.PDE.DNDy(w(:,k,si),this.y(:,k,si-1),tkc,yk,1);
                                DtDmtw = DtDmtw - this.DexpDm_times(this.c(si),v,yk,dtk,yk,1);
                                for sj = 1:si-1
                                    DtDmtw = DtDmtw - dtk*this.DaDm_times(si,sj,v,this.Y(:,k,sj),dtk,yk,1);
                                end
                            end
                        end
                    end
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
        
        % ==========================================================
        
        function results = run_tests(this,varargin)
            
            test_phi_derivatives = 0;
            test_phi_transpose   = 0;
            
            results = [];

            this.precompute_coefficients();
            
            dt = 0.25;

            N  = this.PDE.get('N');
            Nm = this.PDE.get('Nm');
            K  = length(this.t)-1;

            q  = 0.05*rand(N,K,this.s+1);
            y  = this.compute_T_inv(q);
            yk = y(:,10,this.s+1);
            
            if test_phi_derivatives

                fprintf('Testing derivatives of phi:');
                
                dyk = 0.1*norm(yk)*rand(N,1);
                w   = rand(N,1);

                phi               = this.phi_times(2,2,w,dt,yk);
                DphiDyk_times_dyk = this.DphiDyk_times(2,2,dyk,w,dt,yk);
                a = zeros(N,1);
                DaDyk_times_dyk = zeros(N,1);
                for si = 1:this.s
                    for sj = 1:si-1
                        a               = a + this.a_times_test(si,sj,w,dt,yk);
                        DaDyk_times_dyk = DaDyk_times_dyk + this.DaDyk_times(si,sj,dyk,w,dt,yk);
                    end
                end
                b                 = this.b_times_test(4,w,dt,yk);
                DbDyk_times_dyk   = this.DbDyk_times(4,dyk,w,dt,yk);
                e                 = this.exp_times_test(1,w,dt,yk);
                DexpDyk_times_dyk = this.DexpDyk_times(1,dyk,w,dt,yk);

                for i = 1:10
                    i
                    h = 1*0.1^(i+1);

                    phi_pert = this.phi_times(2,2,w,dt,yk+h*dyk);
                    a_pert = zeros(N,1);
                    for si = 1:this.s
                        for sj = 1:si-1
                            a_pert = a_pert + this.a_times_test(si,sj,w,dt,yk+h*dyk);
                        end
                    end
                    b_pert   = this.b_times_test(4,w,dt,yk+h*dyk);
                    exp_pert = this.exp_times_test(1,w,dt,yk+h*dyk);

                    err_dphi(i)    = norm(phi_pert(:)-phi(:))/norm(phi(:));
                    err_dphidyk(i) = norm(phi_pert(:)-phi(:)-h*DphiDyk_times_dyk(:))/norm(phi(:));
                    err_da(i)      = norm(a_pert(:)-a(:))/norm(a(:));
                    err_dadyk(i)   = norm(a_pert(:)-a(:)-h*DaDyk_times_dyk(:))/norm(a(:));
                    err_db(i)      = norm(b_pert(:)-b(:))/norm(b(:));
                    err_dbdyk(i)   = norm(b_pert(:)-b(:)-h*DbDyk_times_dyk(:))/norm(b(:));
                    err_dexp(i)    = norm(exp_pert(:)-e(:))/norm(e(:));
                    err_dexpdyk(i) = norm(exp_pert(:)-e(:)-h*DexpDyk_times_dyk(:))/norm(e(:));
                end
 
                results.phi.dphi    = err_dphi';
                results.phi.dphidyk = err_dphidyk';
                results.a.da        = err_da';
                results.a.dadyk     = err_dadyk';
                results.b.db        = err_db';
                results.b.dbdyk     = err_dbdyk';
                results.exp.dexp    = err_dexp';
                results.exp.dexpdyk = err_dexpdyk';
            elseif test_phi_transpose
                fprintf('Testing transposes of phi:');
                
                err_DtDm_t = 0;
                err_DtDyk_t = 0;
                
                w = rand(Nm,1);
                z = rand(N,1);
                u = rand(N,1);
                v = rand(N,1);
      
                for si = 1:this.s
                    vt = this.DphiDm_times(si,0.5,w,u,dt,yk,0);
                    wt = this.DphiDm_times(si,0.5,v,u,dt,yk,1);

                    err_DtDm_t = err_DtDm_t + (vt(:)'*v(:) - wt(:)'*w(:));
                    
                    vt = this.DphiDm_times(si,1,w,u,dt,yk,0);
                    wt = this.DphiDm_times(si,1,v,u,dt,yk,1);

                    err_DtDm_t = err_DtDm_t + (vt(:)'*v(:) - wt(:)'*w(:));
                end
                
                for si = 1:this.s
                    vt = this.DphiDyk_times(si,0.5,z,u,dt,yk,0);
                    zt = this.DphiDyk_times(si,0.5,v,u,dt,yk,1);

                    err_DtDyk_t = err_DtDyk_t + (vt(:)'*v(:) - zt(:)'*z(:))/(vt(:)'*v(:));
                    
                    vt = this.DphiDyk_times(si,1,z,u,dt,yk,0);
                    zt = this.DphiDyk_times(si,1,v,u,dt,yk,1);
                    
                    err_DtDyk_t = err_DtDyk_t + (vt(:)'*v(:) - zt(:)'*z(:))/(vt(:)'*v(:));
                end
                
                err_DtDm_t
                err_DtDyk_t
            else
                results = run_tests@Integrator(this,varargin{:});
            end
        end
    end
    
    methods (Access = private)
        
        function precompute_coefficients(this)
            
            M  = this.contour.M;
            z  = this.contour.z(:);
            zp = this.contour.zp(:);

            exp_z = exp(z);

            this.a_coeffs = zeros(M,this.s,this.s,2);
            this.b_coeffs = zeros(M,this.s);

            if this.scheme == 1              % Euler
                % Set up b_coeffs
                this.b_coeffs(:,1) = (exp_z./z).*zp;
            elseif this.scheme == 2          % Cox-Matthews
                % Set up a_coeffs
                this.a_coeffs(:,2,1,1) = (0.5*exp_z./z).*zp;
                this.a_coeffs(:,3,2,1) = (0.5*exp_z./z).*zp;
                this.a_coeffs(:,4,1,1) = (-exp_z./z).*zp;
                this.a_coeffs(:,4,1,2) = (exp_z./z).*zp;
                this.a_coeffs(:,4,3,1) = (exp_z./z).*zp;

                % Set up b_coeffs
                this.b_coeffs(:,1) = (exp_z./z - 3*exp_z./z.^2 + 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,2) = (2*exp_z./z.^2 - 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,3) = (2*exp_z./z.^2 - 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,4) = (4*exp_z./z.^3 - exp_z./z.^2).*zp;
            elseif this.scheme == 3          % Krogstad
                % Set up a_coeffs
                this.a_coeffs(:,2,1,1) = (0.5*exp_z./z).*zp;
                this.a_coeffs(:,3,1,1) = (0.5*exp_z./z - exp_z./z.^2).*zp;
                this.a_coeffs(:,3,2,1) = (exp_z./z.^2).*zp;
                this.a_coeffs(:,4,1,2) = (exp_z./z - 2*exp_z./z.^2).*zp;
                this.a_coeffs(:,4,3,2) = (2*exp_z./z.^2).*zp;

                % Set up b_coeffs
                this.b_coeffs(:,1) = (exp_z./z - 3*exp_z./z.^2 + 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,2) = (2*exp_z./z.^2 - 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,3) = (2*exp_z./z.^2 - 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,4) = (4*exp_z./z.^3 - exp_z./z.^2).*zp;
            elseif this.scheme == 4          % Hochbruck-Ostermann
                % Set up a_coeffs
                this.a_coeffs(:,2,1,1) = (0.5*exp_z./z).*zp;
                this.a_coeffs(:,3,1,1) = (0.5*exp_z./z - exp_z./z.^2).*zp;
                this.a_coeffs(:,3,2,1) = (exp_z./z.^2).*zp;
                this.a_coeffs(:,4,1,2) = (exp_z./z - 2*exp_z./z.^2).*zp;
                this.a_coeffs(:,4,2,2) = (exp_z./z.^2).*zp;
                this.a_coeffs(:,4,3,2) = (exp_z./z.^2).*zp;
                this.a_coeffs(:,5,1,1) = (0.5*exp_z./z - 0.75*exp_z./z.^2 + 0.5*exp_z./z.^3).*zp;
                this.a_coeffs(:,5,1,2) = (-0.25*exp_z./z.^2 + exp_z./z.^3).*zp;
                this.a_coeffs(:,5,2,1) = (0.5*exp_z./z.^2 - 0.5*exp_z./z.^3).*zp;
                this.a_coeffs(:,5,2,2) = (0.25*exp_z./z.^2 - exp_z./z.^3).*zp;
                this.a_coeffs(:,5,3,1) = (0.5*exp_z./z.^2 - 0.5*exp_z./z.^3).*zp;
                this.a_coeffs(:,5,3,2) = (0.25*exp_z./z.^2 - exp_z./z.^3).*zp;
                this.a_coeffs(:,5,4,1) = (0.5*exp_z./z.^3 - 0.25*exp_z./z.^2).*zp;
                this.a_coeffs(:,5,4,2) = (-0.25*exp_z./z.^2 + exp_z./z.^3).*zp;

                % Set up b_coeffs
                this.b_coeffs(:,1) = (exp_z./z - 3*exp_z./z.^2 + 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,4) = (-exp_z./z.^2 + 4*exp_z./z.^3).*zp;
                this.b_coeffs(:,5) = (4*exp_z./z.^2 - 8*exp_z./z.^3).*zp;
            end

            this.a_coeffs = (-1i/M)*this.a_coeffs;
            this.b_coeffs = (-1i/M)*this.b_coeffs;

            % Set up exp_coeffs
            this.exp_coeffs = (-1i/M)*exp_z.*zp;
        end
        
        function phiijw = phi_times(this,i,j,w,dt,varargin)
            
            M = this.contour.M;
            z = this.contour.z;
            N = this.PDE.get('N');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            phiijw = zeros(N,1);
            
            if this.PDE.L_diagonal()
                if isempty(this.R05) || isempty(this.R1)
                    this.prepare_resolvents(dt);
                end
                for m = 1:M
                    v05(:,m) = w.*this.R05(:,m);
                    v1(:,m)  = w.*this.R1(:,m);
                    
                    if j == 2 || j == 3
                        phiijw = phiijw + (exp(z(m))/z(m)^i)*v05(:,m);
                    elseif j == 4
                        phiijw = phiijw + (exp(z(m))/z(m)^i)*v1(:,m);
                    end
                end
            else
                yk = []; tp = 0;
                if nargin > 5
                    yk = varargin{1};
                    if nargin > 6, tp = varargin{2}; end
                end

                for m = 1:M
                    R05 = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s -     dt*this.PDE.L(s,yk,tp));
                    
                    v05(:,m) = this.solve_resolvent(R05,w);
                    v1(:,m)  = this.solve_resolvent(R1,w);
                    
                    if j == 2 || j == 3
                        phiijw = phiijw + (exp(z(m))/z(m)^i)*v05(:,m);
                    elseif j == 4
                        phiijw = phiijw + (exp(z(m))/z(m)^i)*v1(:,m);
                    end
                end
            end
            
            phiijw = (1/M)*phiijw;
            if ~this.PDE.L_diagonal()
                phiijw = real(phiijw);
            end
        end
        
        function aijw = a_times(this,i,j)
            aijw = (this.u05*this.a_coeffs(:,i,j,1) + this.u1*this.a_coeffs(:,i,j,2));
            if ~this.PDE.L_diagonal()
                aijw = real(aijw);
            end
        end
        
        function aijw = a_times_test(this,i,j,w,dt,varargin)
            
            M = this.contour.M;
            z = this.contour.z;
            N = this.PDE.get('N');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            aijw = zeros(N,1);

            yk = []; tp = 0;
            if nargin > 5
                yk = varargin{1};
                if nargin > 6, tp = varargin{2}; end
            end
            
            for m = 1:M
                R05 = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,tp));
                R1  = @(s)(z(m)*s -     dt*this.PDE.L(s,yk,tp));

                v05(:,m) = this.solve_resolvent(R05,w);
                v1(:,m)  = this.solve_resolvent(R1,w);

                aijw = aijw + this.a_coeffs(m,i,j,1)*v05(:,m);
                aijw = aijw + this.a_coeffs(m,i,j,2)*v1(:,m);
            end
            
            if ~this.PDE.L_diagonal()
                aijw = real(aijw);
            end
        end
        
        function biw = b_times(this,i)
            biw = (this.u1*this.b_coeffs(:,i));
            if ~this.PDE.L_diagonal()
                biw = real(biw);
            end
        end
        
        function biw = b_times_test(this,i,w,dt,varargin)
            
            M = this.contour.M;
            z = this.contour.z;
            N = this.PDE.get('N');
            
            v   = zeros(N,M);
            biw = zeros(N,1);

            yk = []; tp = 0;
            if nargin > 4
                yk = varargin{1};
                if nargin > 5, tp = varargin{2}; end
            end
            
            for m = 1:M
                R = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));
                v(:,m) = this.solve_resolvent(R,w);
                biw = biw + this.b_coeffs(m,i)*v(:,m);
            end
            
            if ~this.PDE.L_diagonal()
                biw = real(biw);
            end
        end
        
        function expw = exp_times(this,c)
            if c == 0.5
                expw = (this.u05*this.exp_coeffs);
            elseif c == 1
                expw = (this.u1*this.exp_coeffs);
            end
            if ~this.PDE.L_diagonal()
                expw = real(expw);
            end
        end
        
        function expw = exp_times_test(this,c,w,dt,varargin)
            
            M = this.contour.M;
            z = this.contour.z;
            N = this.PDE.get('N');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 4
                yk = varargin{1};
                if nargin > 5, tp = varargin{2}; end
            end

            for m = 1:M
                R05 = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,tp));
                R1  = @(s)(z(m)*s -     dt*this.PDE.L(s,yk,tp));

                v05(:,m) = this.solve_resolvent(R05,w);
                v1(:,m)  = this.solve_resolvent(R1,w);
            end
            
            if c == 0.5
                expw = (v05*this.exp_coeffs);
            elseif c == 1
                expw = (v1*this.exp_coeffs);
            end
            if ~this.PDE.L_diagonal()
                expw = real(expw);
            end
        end

        function aijw = DaDyk_times(this,i,j,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            DvDyk05 = zeros(N,M);
            DvDyk1  = zeros(N,M);
                
            yk = []; tp = 0;
            if nargin > 6
                yk = varargin{1};
                if nargin > 7, tp = varargin{2}; end
            end

            if tp == 1
                
                for m = 1:M
                    R05          = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk));
                    v05(:,m)     = this.solve_resolvent(R05,w);
                    R05          = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,1));
                    temp         = this.solve_resolvent(R05,u);
                    DvDyk05(:,m) = 0.5*dt*this.PDE.DLDyk(temp,v05(:,m),yk,tp);
                    
                    R1           = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v1(:,m)      = this.solve_resolvent(R1,w);
                    R1           = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp         = this.solve_resolvent(R1,u);
                    DvDyk1(:,m)  = dt*this.PDE.DLDyk(temp, v1(:,m),yk,tp);
                end
            else
                DvDyk05 = zeros(N,M);
                DvDyk1  = zeros(N,M);
                
                for m = 1:M
                    R05 = @(s)(z(m)*s - (dt/2)*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));

                    v05(:,m) = this.solve_resolvent(R05,w);
                    v1(:,m)  = this.solve_resolvent(R1,w);

                    DvDyk05(:,m) = this.solve_resolvent(R05,0.5*dt*this.PDE.DLDyk(u,v05(:,m),yk,tp));
                    DvDyk1(:,m)  = this.solve_resolvent(R1,     dt*this.PDE.DLDyk(u,v1(:,m),yk,tp));
                end
            end
            
            aijw = (DvDyk05*this.a_coeffs(:,i,j,1) + DvDyk1*this.a_coeffs(:,i,j,2));
            if ~this.PDE.L_diagonal()
                aijw = real(aijw);
            end
        end
        
        function biw = DbDyk_times(this,i,u,w,dt,varargin)
            
            M = this.contour.M;
            z = this.contour.z;
            N = this.PDE.get('N');
            
            v     = zeros(N,M);
            DvDyk = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 5
                yk = varargin{1};
                if nargin > 6, tp = varargin{2}; end
            end

            if tp == 1
                for m = 1:M
                    R          = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v(:,m)     = this.solve_resolvent(R,w);
                    R          = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp       = this.solve_resolvent(R,u);
                    DvDyk(:,m) = dt*this.PDE.DLDyk(temp,v(:,m),yk,tp);
                end
            else
                for m = 1:M
                    R          = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));
                    v(:,m)     = this.solve_resolvent(R,w);
                    DvDyk(:,m) = this.solve_resolvent(R,dt*this.PDE.DLDyk(u,v(:,m),yk,tp));
                end
            end
            
            biw = (DvDyk*this.b_coeffs(:,i));
            if ~this.PDE.L_diagonal()
                biw = real(biw);
            end
        end
        
        function expw = DexpDyk_times(this,c,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            DvDyk05 = zeros(N,M);
            DvDyk1  = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 5
                yk = varargin{1};
                if nargin > 6, tp = varargin{2}; end
            end

            if tp == 1
                for m = 1:M
                    R05          = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk));
                    v05(:,m)     = this.solve_resolvent(R05,w);
                    R05          = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,1));
                    temp         = this.solve_resolvent(R05,u);
                    DvDyk05(:,m) = 0.5*dt*this.PDE.DLDyk(temp,v05(:,m),yk,tp);
                    
                    R1           = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v1(:,m)      = this.solve_resolvent(R1,w);
                    R1           = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp         = this.solve_resolvent(R1,u);
                    DvDyk1(:,m)  = dt*this.PDE.DLDyk(temp, v1(:,m),yk,tp);
                end
            else
                for m = 1:M
                    R05 = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s -     dt*this.PDE.L(s,yk,tp));

                    v05(:,m) = this.solve_resolvent(R05,w);
                    v1(:,m)  = this.solve_resolvent(R1,w);

                    DvDyk05(:,m) = this.solve_resolvent(R05,0.5*dt*this.PDE.DLDyk(u,v05(:,m),yk,tp));
                    DvDyk1(:,m)  = this.solve_resolvent(R1,     dt*this.PDE.DLDyk(u, v1(:,m),yk,tp));
                end
            end
            
            if c == 0.5
                expw = (DvDyk05*this.exp_coeffs);
            elseif c == 1
                expw = (DvDyk1*this.exp_coeffs);
            end
            if ~this.PDE.L_diagonal()
                expw = real(expw);
            end
        end
        
        function phiijw = DphiDyk_times(this,i,c,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            
            v      = zeros(N,M);
            DvDyk  = zeros(N,M);
            phiijw = zeros(N,1);
            
            yk = []; tp = 0;
            if nargin > 6
                yk = varargin{1};
                if nargin > 7, tp = varargin{2}; end
            end

            if tp == 1
                for m = 1:M
                    R          = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk));
                    v(:,m)     = this.solve_resolvent(R,w);
                    R          = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk,1));
                    temp       = this.solve_resolvent(R,u);
                    DvDyk(:,m) = c*dt*this.PDE.DLDyk(temp,v(:,m),yk,1);
                end
            else
                for m = 1:M
                    R = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk,tp));
                    v(:,m) = this.solve_resolvent(R,w);
                    DvDyk(:,m) = this.solve_resolvent(R,c*dt*this.PDE.DLDyk(u,v(:,m),yk,tp));
                end
            end
            
            for m = 1:M
                phiijw = phiijw + (exp(z(m))/z(m)^i)*DvDyk(:,m);
            end
            
            phiijw = (1/M)*phiijw;
            if ~this.PDE.L_diagonal()
                phiijw = real(phiijw);
            end
        end
        
        function aijw = DaDm_times(this,i,j,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            Nm = this.PDE.get('Nm');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 6
                yk = varargin{1};
                if nargin > 7, tp = varargin{2}; end
            end

            if tp == 1
                DvDm05 = zeros(Nm,M);
                DvDm1  = zeros(Nm,M);
                
                for m = 1:M
                    R05         = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk));
                    v05(:,m)    = this.solve_resolvent(R05,w);
                    R05         = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,1));
                    temp        = this.solve_resolvent(R05,u);
                    DvDm05(:,m) = 0.5*dt*this.PDE.DLDm(temp,v05(:,m),yk,tp);
                    
                    R1          = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v1(:,m)     = this.solve_resolvent(R1,w);
                    R1          = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp        = this.solve_resolvent(R1,u);
                    DvDm1(:,m)  = dt*this.PDE.DLDm(temp, v1(:,m),yk,tp);
                end
            else
                DvDm05 = zeros(N,M);
                DvDm1  = zeros(N,M);
                
                for m = 1:M
                    R05 = @(s)(z(m)*s - (dt/2)*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));

                    v05(:,m) = this.solve_resolvent(R05,w);
                    v1(:,m)  = this.solve_resolvent(R1,w);

                    DvDm05(:,m) = this.solve_resolvent(R05,0.5*dt*this.PDE.DLDm(u,v05(:,m),yk,tp));
                    DvDm1(:,m)  = this.solve_resolvent(R1,     dt*this.PDE.DLDm(u,v1(:,m),yk,tp));
                end
            end
            
            aijw = (DvDm05*this.a_coeffs(:,i,j,1) + DvDm1*this.a_coeffs(:,i,j,2));
            if ~this.PDE.L_diagonal()
                aijw = real(aijw);
            end
        end
        
        function biw = DbDm_times(this,i,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            
            v    = zeros(this.PDE.get('N'),M);
            
            yk = []; tp = 0;
            if nargin > 5
                yk = varargin{1};
                if nargin > 6, tp = varargin{2}; end
            end

            if tp == 1
                DvDm = zeros(this.PDE.get('Nm'),M);
                for m = 1:M
                    R         = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v(:,m)    = this.solve_resolvent(R,w);
                    R         = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp      = this.solve_resolvent(R,u);
                    DvDm(:,m) = dt*this.PDE.DLDm(temp,v(:,m),yk,tp);
                end
            else
                DvDm = zeros(this.PDE.get('N'),M);
                for m = 1:M
                    R         = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));
                    v(:,m)    = this.solve_resolvent(R,w);
                    DvDm(:,m) = this.solve_resolvent(R,dt*this.PDE.DLDm(u,v(:,m),yk,tp));
                end
            end
            
            biw = (DvDm*this.b_coeffs(:,i));
            if ~this.PDE.L_diagonal()
                biw = real(biw);
            end
        end
        
        function expw = DexpDm_times(this,c,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            Nm = this.PDE.get('Nm');
            
            v05 = zeros(N,M);
            v1  = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 5
                yk = varargin{1};
                if nargin > 6, tp = varargin{2}; end
            end

            if tp == 1
                DvDm05 = zeros(Nm,M);
                DvDm1  = zeros(Nm,M);
                
                for m = 1:M
                    R05         = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk));
                    v05(:,m)    = this.solve_resolvent(R05,w);
                    R05         = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,1));
                    temp        = this.solve_resolvent(R05,u);
                    DvDm05(:,m) = 0.5*dt*this.PDE.DLDm(temp,v05(:,m),yk,tp);
                    
                    R1          = @(s)(z(m)*s - dt*this.PDE.L(s,yk));
                    v1(:,m)     = this.solve_resolvent(R1,w);
                    R1          = @(s)(z(m)*s - dt*this.PDE.L(s,yk,1));
                    temp        = this.solve_resolvent(R1,u);
                    DvDm1(:,m)  = dt*this.PDE.DLDm(temp, v1(:,m),yk,tp);
                end
            else
                DvDm05 = zeros(N,M);
                DvDm1  = zeros(N,M);
                
                for m = 1:M
                    R05 = @(s)(z(m)*s - 0.5*dt*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s -     dt*this.PDE.L(s,yk,tp));

                    v05(:,m) = this.solve_resolvent(R05,w);
                    v1(:,m)  = this.solve_resolvent(R1,w);

                    DvDm05(:,m) = this.solve_resolvent(R05,0.5*dt*this.PDE.DLDm(u,v05(:,m),yk,tp));
                    DvDm1(:,m)  = this.solve_resolvent(R1,     dt*this.PDE.DLDm(u, v1(:,m),yk,tp));
                end
            end
            
            if c == 0.5
                expw = (DvDm05*this.exp_coeffs);
            elseif c == 1
                expw = (DvDm1*this.exp_coeffs);
            end
            if ~this.PDE.L_diagonal()
                expw = real(expw);
            end
        end
        
        function phiijw = DphiDm_times(this,i,c,u,w,dt,varargin)
            
            M  = this.contour.M;
            z  = this.contour.z;
            N  = this.PDE.get('N');
            Nm = this.PDE.get('Nm');
            
            v = zeros(N,M);
            
            yk = []; tp = 0;
            if nargin > 6
                yk = varargin{1};
                if nargin > 7, tp = varargin{2}; end
            end

            if tp == 1
                DvDm = zeros(Nm,M);
                phiijw = zeros(Nm,1);
                
                for m = 1:M
                    R         = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk));
                    v(:,m)    = this.solve_resolvent(R,w);
                    R         = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk,1));
                    temp      = this.solve_resolvent(R,u);
                    DvDm(:,m) = c*dt*this.PDE.DLDm(temp,v(:,m),yk,1);
                end
            else
                DvDm = zeros(N,M);
                phiijw = zeros(N,1);
                
                for m = 1:M
                    R         = @(s)(z(m)*s - c*dt*this.PDE.L(s,yk,tp));
                    v(:,m)    = this.solve_resolvent(R,w);
                    DvDm(:,m) = this.solve_resolvent(R,c*dt*this.PDE.DLDm(u,v(:,m),yk,tp));
                end
            end
            
            for m = 1:M
                phiijw = phiijw + (exp(z(m))/z(m)^i)*DvDm(:,m);
            end
            
            phiijw = (1/M)*phiijw;
            if ~this.PDE.L_diagonal()
                phiijw = real(phiijw);
            end
        end
        
        % -----------------------------------------------------------------

        function prepare_resolvents(this,dt,varargin)

            if this.PDE.L_diagonal()
                
                M  = this.contour.M;
                z  = this.contour.z;
                
                if isempty(this.R05) || isempty(this.R1)
                    this.R05 = zeros(this.PDE.get('N'),M);
                    this.R1  = zeros(this.PDE.get('N'),M);
                end
                
                dtL = dt*this.PDE.get_L();
                for m = 1:M
                    this.R05(:,m) = 1./(z(m) - 0.5*dtL);
                    this.R1(:,m)  = 1./(z(m) - dtL);
                end
            else
                this.R05 = [];
                this.R1  = [];
            end
        end
        
        function v = solve_resolvents(this,v,dt,varargin)
            
            this.i2 = this.i2 + 1;
            
            M  = this.contour.M;
            z  = this.contour.z;
            
            if isempty(this.u05) || isempty(this.u1)
                this.u05 = zeros(this.PDE.get('N'),M);
                this.u1  = zeros(this.PDE.get('N'),M);
            end
            
            if this.PDE.L_diagonal()
                if isempty(this.R05) || isempty(this.R1)
                    this.prepare_resolvents(dt);
                end
                for m = 1:M
                    this.u05(:,m) = v.*this.R05(:,m);
                    this.u1(:,m)  = v.*this.R1(:,m);
                end
            else
                yk = []; tp = 0;
                if nargin > 3
                    yk = varargin{1};
                    if nargin > 4, tp = varargin{2}; end
                end

                for m = 1:M
                    R05 = @(s)(z(m)*s - (dt/2)*this.PDE.L(s,yk,tp));
                    R1  = @(s)(z(m)*s - dt*this.PDE.L(s,yk,tp));
                    
                    this.u05(:,m) = this.solve_resolvent(R05,v);
                    this.u1(:,m)  = this.solve_resolvent(R1,v);
                end
            end
        end
        
        function v = solve_resolvent(this,R,w)

            [v,flag,relres,iter] = minres(R,w,1e-8,200);
        end
        
        function setup_contour(this,M,shape)

            if strcmpi(shape,'parabolic')
                theta = pi*(1:2:M-1)/M;
                theta = [-flipud(theta'); theta'];
                z     = M*(.1309-.1194*theta.^2+.2500i*theta);  % Quadrature points along contour
                zp    = M*(-.1194*2*theta+.2500i);              % Derivatives
            elseif strcmpi(shape,'hyperbolic')
                theta = pi*(1:2:M-1)/M;
                theta = [-flipud(theta'); theta'];
                z     = 2.246*M*(1-sin(1.17212 - 0.3443i*theta));  % Quadrature points along contour
                zp    = 0.7732978i*M*(cos(1.17212 - 0.3443i*theta));  % Derivatives
            elseif strcmpi(shape,'circle')
                z     = sr*exp(1i*pi*((1:M)-.5)/M);      % Quadrature points along complex circle
            end
            
            this.contour.M  = M;
            this.contour.z  = z;
            this.contour.zp = zp;
        end
    end
end
