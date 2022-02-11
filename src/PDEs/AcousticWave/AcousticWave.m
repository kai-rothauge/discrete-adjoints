classdef AcousticWave < PDE
    
    properties (SetAccess = private)
        
        sd = 0;             % Spatial discretization:
                            % 2 - 2nd
                            % 4 - 4th
                            % 6 - 6th
                            
        Fourier_basis = 0;  % Only if using a pseudospectral method with
                            % the non-rosenbrock version of the equation
    end
    
    methods (Access = public)
        function this = AcousticWave(varargin)
            this.name      = 'Acoustic Wave';
            this.linearity = 'linear';
            
            s = settings('PDE');
            this.sd         = s.sd;
            this.set_basis();
            
            if nargin == 1 && isnumeric(varargin{1})
                this.set_grid(settings('grid',varargin{1}));
            else
                this.set_grid(settings('grid'));
            end
  
            if nargin > 1
                this.set(varargin{:});
            end
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this, varargin)
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('Options for ''set'':');
                    disp('  ''source'' or ''source term'' or ''q''');
                    disp('  ''initial condition'' or ''IC''');
                    disp('  ''parameters'' or ''m''');
                    disp('  ''grid''');
                    disp('  ''spatial discretization'' or ''sd''');
                elseif strcmpi(parameter,'source') || strcmpi(parameter,'source term') || strcmpi(parameter,'q')
                    this.set_q(varargin{i+1});
                elseif strcmpi(parameter,'initial condition') || strcmpi(parameter,'IC')
                    this.set_IC(varargin{i+1});
                elseif strcmpi(parameter,'parameters') || strcmpi(parameter,'m')
                    this.set_parameters(varargin{i+1});
                elseif strcmpi(parameter,'grid')
                    this.set_grid(varargin{i+1});
                elseif strcmpi(parameter,'spatial discretization') || strcmpi(parameter,'sd')
                    this.set_sd(varargin{i+1});
                end
            end
        end
        
        function set_source(this,q)
            this.q = q;
        end
        
        function set_q(this,q)
            this.q = q;
        end
        
        function set_IC(this,IC)
            if ischar(IC)
                this.setup_IC(IC);
            else
                this.u0 = IC;
                if this.Fourier_basis
                    this.u0 = this.change_basis(this.u0);
                end
            end
        end
        
        function set_parameters(this,m)
            if ischar(m)
                this.setup_parameters(m);
            else
                sets = this.split_parameters(m);

                this.m.r = sets{1};
                this.m.g = sets{2};
            end
        end

        function set_grid(this,grid)
            this.grid = grid;
            
            this.setup_IC();
            this.setup_parameters();
        end
        
        function set_sd(this,sd)
            if sd == 6 || strcmpi(sd,'6') || strcmpi(sd,'6th order') || strcmpi(sd,'6th')
                this.sd = 6;
            elseif sd == 4 || strcmpi(sd,'4') || strcmpi(sd,'4th order') || strcmpi(sd,'4th')
                this.sd = 4;
            else
                this.sd = 2;
            end
            
            this.setup_Lhat();
            this.set_basis();
        end

        % ---------------------- Accessor functions -----------------------
        
        function out = get(this,parameter)
            out = [];
            
            if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                disp('Options for ''get'':');
                disp('  ''num unknowns'' or ''N'' or ''number of unknowns''');
                disp('  ''num parameters'' or ''Nm'' or ''number of parameters''');
                disp('  ''initial condition'' or ''IC''');
                disp('  ''parameters'' or ''m''');
                disp('  ''spatial discretization'' or ''sd''');
                disp('  ''grid''');
            elseif strcmpi(parameter,'num unknowns') || strcmpi(parameter,'N') ...
                                        || strcmpi(parameter,'number of unknowns')
                out = this.grid.x.N*this.grid.y.N;
            elseif strcmpi(parameter,'num parameters') || strcmpi(parameter,'Nm') ...
                                        || strcmpi(parameter,'number of parameters')
                out = this.Nm.r + this.Nm.g;
            elseif strcmpi(parameter,'initial condition') || strcmpi(parameter,'IC')
                out = this.u0;
            elseif strcmpi(parameter,'parameters') || strcmpi(parameter,'m')
                out = [this.m.r; this.m.g];
            elseif strcmpi(parameter,'random parameters') || strcmpi(parameter,'dm')
                out = [norm(this.m.r(:))*rand(length(this.m.r(:)),1);
                     norm(this.m.g(:))*rand(length(this.m.g(:)),1)];
            elseif strcmpi(parameter,'domain dimensions')
                out = [this.grid.x.N; this.grid.y.N];
            elseif strcmpi(parameter,'spatial discretization') || strcmpi(parameter,'sd')
                out = this.sd;
            elseif strcmpi(parameter,'grid')
                out = this.grid;
            end
        end

        function m = get_parameters(this)
            m = this.m;
        end
        
        function grid = get_grid(this)
            grid = this.grid;
        end
        
        function u0 = get_IC(this)
            u0 = this.u0;
%             if this.Fourier_basis
%                 u0 = this.change_basis(u0,'inverse');
%             end
        end
        
        function N = get_N(this)
            N = this.grid.x.N*this.grid.y.N;
        end
        
        function Nm = get_num_parameters(this)
            Nm = this.Nm.r + this.Nm.g;
        end
        
        function m = get_parameter_names(this)
            m{1}.name = 'r';
            m{2}.name = 'g';
        end
        
        function Nm = get_Nm(this)
            Nm = this.Nm.r + this.Nm.g;
        end
        
        function m_grad = get_parameter_gradient(this)
            
            Nx = this.grid.x.N;
            Ny = this.grid.y.N;
            
            xx  = linspace(this.grid.x.min,this.grid.x.max,Nx+1); 
            yy  = linspace(this.grid.y.min,this.grid.y.max,Ny+1);
            
            dx = xx(2)-xx(1);
            dy = yy(2)-yy(1);

            D2x = this.derivative_1D(2,Nx,dx);
            D2y = this.derivative_1D(2,Ny,dy);

            D2 = kron(D2x,speye(Ny)) + kron(speye(Nx),D2y);
            
            m_grad = [D2 sparse(Nx*Ny,Nx*Ny);
                      sparse(Nx*Ny,Nx*Ny) D2];
        end
        
        % -----------------------------------------------------------------

        function display_parameters(this,varargin)
            
            if nargin > 2
                m_sets = this.split_parameters(varargin{1});
                m.r = m_sets{1};
                m.g = m_sets{2};
                ah = varargin{2};
            else
                m = this.m;
                ah = varargin{1};
            end
            
            if ~isempty(ah{1}) || ~isempty(ah{2})
            
                x = linspace(this.grid.x.min,this.grid.x.max,this.grid.x.N+1); x = x(1:end-1);
                y = linspace(this.grid.y.min,this.grid.y.max,this.grid.y.N+1); y = y(1:end-1);

                if ~isempty(ah{1})
                    imagesc(x,y,reshape(m.r,this.grid.x.N,this.grid.y.N),'Parent',ah{1});
                    colormap(ah{1},'parula')
                    colorbar(ah{1})
                    xlabel(ah{1},'x','Fontsize',10)
                    ylabel(ah{1},'y','rotation',0,'Fontsize',10)
                    axis(ah{1},'equal','tight')
                    set(ah{1},'CLim',[this.m_bounds.r.min,this.m_bounds.r.max])
                    drawnow
                end

                if ~isempty(ah{2})
                    imagesc(x,y,reshape(m.g,this.grid.x.N,this.grid.y.N),'Parent',ah{2});
                    colormap(ah{2},'parula')
                    colorbar(ah{2})
                    xlabel(ah{2},'x','Fontsize',10)
                    ylabel(ah{2},'y','rotation',0,'Fontsize',10)
                    axis(ah{2},'equal','tight')
                    set(ah{2},'CLim',[this.m_bounds.g.min,this.m_bounds.g.max])
                    drawnow
                end
            end
        end
        
        function display_initial_condition(this,ah)
            
            if ~isempty(ah)
                if this.Fourier_basis
                    u0 = this.change_basis(this.u0,'inverse');
                else
                    u0 = this.u0;
                end
            
                x = linspace(this.grid.x.min,this.grid.x.max,this.grid.x.N+1); x = x(1:end-1);
                y = linspace(this.grid.y.min,this.grid.y.max,this.grid.y.N+1); y = y(1:end-1);

                imagesc(x,y,reshape(u0,this.grid.x.N,this.grid.y.N),'Parent',ah);
                colormap(ah,'parula')
                colorbar(ah)
                title(ah,'Initial condition','Fontsize',10);
                xlabel(ah,'x','Fontsize',10)
                ylabel(ah,'y','rotation',0,'Fontsize',10)
                axis(ah,'equal','tight')
                drawnow
            end
        end
        
        function display_solution(this,u,t,ah)
            
            if ~isempty(ah)
                if this.Fourier_basis
                    u = this.change_basis(u,'inverse');
                end

                x = linspace(this.grid.x.min,this.grid.x.max,this.grid.x.N+1); x = x(1:end-1);
                y = linspace(this.grid.y.min,this.grid.y.max,this.grid.y.N+1); y = y(1:end-1);

                imagesc(x,y,reshape(u,this.grid.x.N,this.grid.y.N),'Parent',ah);
                colormap(ah,'parula')
                colorbar(ah)
                str = sprintf('Solution of %s equation at time t = %5.2fs',this.name,t);
                title(ah,str,'Fontsize',10);
                xlabel(ah,'x','Fontsize',10)
                ylabel(ah,'y','rotation',0,'Fontsize',10)
                axis(ah,'equal','tight')
                drawnow
            end
        end
        
        function sets = split_parameters(this,varargin)
            
            if nargin > 1
                w_m = varargin{1};
                sets{1} = w_m(1:this.Nm.r);
                sets{2} = w_m(this.Nm.r+1:end);
            else
                sets{1} = m.r; 
                sets{2} = m.g;
            end
        end
        
        function m = update_parameters(this,dm)
            
            dm = this.split_parameters(dm);
            
            this.m.r = this.m.r + dm{1};
            this.m.g = this.m.g + dm{2};
            
            m = [this.m.r; this.m.g];
        end
        
        function m = constrain_parameters(this,varargin)
            
            if nargin > 1
                m_sets = this.split_parameters(varargin{1});
                m.r = m_sets{1};
                m.g = m_sets{2};
                
                m.r(m.r < this.m_bounds.r.min) = this.m_bounds.r.min;
                m.r(m.r > this.m_bounds.r.max) = this.m_bounds.r.max;
                m.g(m.g < this.m_bounds.g.min) = this.m_bounds.g.min;
                m.g(m.g > this.m_bounds.g.max) = this.m_bounds.g.max;
                
                m = [m.r; m.g];
            else
                this.m.r(this.m.r < this.m_bounds.r.min) = this.m_bounds.r.min;
                this.m.r(this.m.r > this.m_bounds.r.max) = this.m_bounds.r.max;
                this.m.g(this.m.g < this.m_bounds.g.min) = this.m_bounds.g.min;
                this.m.g(this.m.g > this.m_bounds.g.max) = this.m_bounds.g.max;
                
                m = [this.m.r; this.m.g];
            end
        end

        function Ly = L(this,y,yk,varargin)

            Ly = this.Lhat_times(y,varargin{:});
            if this.rosenbrock && ~isempty(yk)
                Ly = Ly + (this.m.r + 2*this.m.g.*yk - 3*(yk.^2)).*y;
            end
        end
        
        function Ny = N(this,y,t,yk)
            
            if this.rosenbrock && ~isempty(yk)
                Ny = this.m.g.*(y.^2 - 2*yk.*y) - (y.^3 - 3*(yk.^2).*y);
            else
                if this.Fourier_basis
                    y = this.change_basis(y,'inverse');
                end
                Ny = this.m.r.*y + this.m.g.*(y.^2) - y.^3;
                if this.Fourier_basis
                    Ny = this.change_basis(Ny);
                end
            end
        end
        
        function diagonal = L_diagonal(this)
            if this.Fourier_basis
                diagonal = 1;
            else
                diagonal = 0;
            end
        end
        
        function u = change_basis(this,u,varargin)

            Nx = this.grid.x.N;
            Ny = this.grid.y.N;
            K = size(u,2);
            
            for k = 1:K
                if nargin > 2
                    utemp = ifft2(reshape(full(u(:,k)),Nx,Ny));
                    utemp = real(utemp);
                else
                    utemp = fft2(reshape(full(u(:,k)),Nx,Ny));
                end
                u(:,k) = (utemp(:));
            end
        end
        
        function L = get_L(this)
            if this.Fourier_basis
                L = this.Lhat(:);
            else
                L = this.Lhat;
            end
        end
        
        % -------------- Derivative of L and N wrt solution ---------------
        
        function dLv = DLDyk(this,v,w,varargin)
            
            yk = []; transpose = 0;
            if nargin > 3
                yk = varargin{1};
                if nargin > 4, transpose = varargin{2}; end
            end
            
            if this.rosenbrock && ~isempty(yk)
                dLv = (2*this.m.g - 6*yk).*v.*w;
            else
                dLv = zeros(length(v),1);
            end
        end
        
        function dN = DNDyk(this,w,v,t,varargin)

            yk = []; transpose = 0;
            if nargin > 4
                yk = varargin{1};
                if nargin > 5, transpose = varargin{2}; end
            end

            dN = [];
            if this.rosenbrock && ~isempty(yk)
                dN = -2*(this.m.g - 3*yk).*v.*w;
            end
        end
        
        function dN = DNDy(this,w,v,t,varargin)

            yk = []; transpose = 0;
            if nargin > 4
                yk = varargin{1};
                if nargin > 5, transpose = varargin{2}; end
            end

            if this.rosenbrock && ~isempty(yk)
                dN = (2*this.m.g.*(v-yk) - 3*(v.^2-yk.^2)).*w;
            else
                if this.Fourier_basis
                    w = this.change_basis(w,'inverse');
                    v = this.change_basis(v,'inverse');
                end
                dN = (this.m.r + 2*this.m.g.*v - 3*(v.^2)).*w;
                if this.Fourier_basis
                    dN = this.change_basis(dN);
                end
            end
        end

        % ------------- Derivative of N wrt model parameters --------------
        
        function dLv = DLDm(this,w,v,varargin)
            
            yk = []; transpose = 0;
            if nargin > 3
                yk = varargin{1};
                if nargin > 4, transpose = varargin{2}; end
            end
            
            if transpose
                if this.rosenbrock && ~isempty(yk)
                    dLv = [v.*w; 2*yk.*v.*w];
                else
                    if this.Fourier_basis
                        v = this.change_basis(v,'inverse');
                        w = numel(w)*this.change_basis(w,'inverse');
                    end
                    dLv = [zeros(size(w)); 2*yk.*v.*w];
                end
            else
                w_sets = this.split_parameters(w);
                w_r = w_sets{1}; w_g = w_sets{2};
                if this.rosenbrock && ~isempty(yk)
                    dLv = (w_r + 2*yk.*w_g).*v;
                else
                    dLv = zeros(length(v),1);
                end
            end
        end
        
        function dN = DNDm(this,w,v,t,varargin)
            
            yk = []; transpose = 0;
            if nargin > 4
                yk = varargin{1};
                if nargin > 5, transpose = varargin{2}; end
            end

            if transpose
                if this.rosenbrock && ~isempty(yk)
                    dN = [zeros(size(w)); (v.^2-2*yk.*v).*w];
                else
                    if this.Fourier_basis
                        v = this.change_basis(v,'inverse');
                        w = numel(w)*this.change_basis(w,'inverse');
                    end
                    dN = [v.*w; v.^2.*w];
                end
            else
                w_sets = this.split_parameters(w);
                w_r = w_sets{1}; w_g = w_sets{2};
                if this.rosenbrock && ~isempty(yk)
                    dN = (v.^2-2*yk.*v).*w_g;
                else
                    if this.Fourier_basis
                        v = this.change_basis(v,'inverse');
                    end
                    dN = v.*w_r + v.^2.*w_g; 
                    if this.Fourier_basis
                        dN = this.change_basis(dN); 
                    end
                end
            end
        end
        
        % ------------------------ Source function ------------------------
        
        function qt = q(this,t)
            qt = 0;
        end
        
        % -------------------------- Diagnostics --------------------------
        
        function run_diagnostics(this,varargin)
            
            if nargin > 1
                selection = varargin{1};
            end
            
            this.run_tests(selection);
        end
        
        function results = run_tests(this)

            disp('Testing N:')

            Nx = this.grid.x.N;
            Ny = this.grid.y.N;
            N  = Nx*Ny;
            
            t   = rand;
            u   = rand(N,1);
            uk  = rand(N,1);
            du  = rand(N,1);
            duk = rand(N,1);
            
            if this.Fourier_basis
                u  = fft2(reshape(u, Nx,Ny)); u  = u(:);
                uk = fft2(reshape(uk,Nx,Ny)); uk = uk(:);
                du = fft2(reshape(du,Nx,Ny)); du = du(:);
            end
            
            old_m = this.m;
            
            dr = rand(size(this.m.r));
            dg = rand(size(this.m.g));
            dm = [dr;dg];

            Nu            = this.N(u,t,uk);
            DNDm_times_dM = this.DNDm(dm,u,t,uk,0);
            DNDy_times_dY = this.DNDy(du,u,t,uk,0);
            
            Lu            = this.L(u,uk);
            DLDm_times_dM = this.DLDm(dm,u,uk,0);
            DLDyk_times_dY = this.DLDyk(duk,u,uk,0);
            
            for i = 1:10
                h = 1*0.1^(i+1);
                
                this.m.r = old_m.r + h*dr;
                this.m.g = old_m.g + h*dg;

                Nu_pert = this.N(u,t,uk);
                Lu_pert = this.L(u,uk);
                
                this.m.r = old_m.r;
                this.m.g = old_m.g;

                err_N(i)    = norm(Nu_pert - Nu)/norm(Nu);
                err_DNDm(i) = norm(Nu_pert - Nu - h*DNDm_times_dM)/norm(Nu);
                err_L(i)    = norm(Lu_pert - Lu)/norm(Lu);
                err_DLDm(i) = norm(Lu_pert - Lu - h*DLDm_times_dM)/norm(Lu);
            end
            
            err_NM   = err_N';
            err_DNDm = err_DNDm';
            err_LM   = err_L';
            err_DLDm = err_DLDm';
            
            for i = 1:10
                h = 1*0.1^(i+1);
                
                Nu_pert = this.N(u+h*du,t,uk);
                Lu_pert = this.L(u,uk+h*duk);

                err_N(i)     = norm(Nu_pert - Nu,2)/norm(Nu,2);
                err_DNDy(i)  = norm(Nu_pert - Nu - h*DNDy_times_dY,2)/norm(Nu,2);
                err_L(i)     = norm(Lu_pert - Lu,2)/norm(Lu,2);
                err_DLDyk(i) = norm(Lu_pert - Lu - h*DLDyk_times_dY,2)/norm(Lu,2);
            end
            
            err_NY    = err_N';
            err_DNDy  = err_DNDy';
            err_LY    = err_L';
            err_DLDyk = err_DLDyk';
            
            results.error_NM   = err_NM;
            results.error_DNDm = err_DNDm;
            results.error_NY   = err_NY;
            results.error_DNDy = err_DNDy;
            results.error_LM   = err_LM;
            results.error_DLDm = err_DLDm;
            results.error_LY   = err_LY;
            results.error_DLDyk = err_DLDyk;

            % ============================ Test transposes ============================

            % Check transpose of Jacobian
            wM = [rand(size(this.m.r));
                  rand(size(this.m.g))];
            wY = rand(N,1);
            
            wtM = rand(N,1);
            wtY = rand(N,1);
            
            if this.Fourier_basis
                wY  = fft2(reshape(wY,Nx,Ny));  wY  = wY(:);
                wtM = fft2(reshape(wtM,Nx,Ny)); wtM = wtM(:);
                wtY = fft2(reshape(wtY,Nx,Ny)); wtY = wtY(:);
            end
            
            vM = this.DNDm(wM,u,t,uk,0);
            vY = this.DNDy(wY,u,t,uk,0);
            
            vtM = this.DNDm(wtM,u,t,uk,1);
            vtY = this.DNDy(wtY,u,t,uk,1);

            err_DNDm_t = ((wtM(:))'*vM(:) - (wM(:))'*vtM(:))/abs((wtM(:))'*vM(:));
            err_DNDy_t = ((wtY(:))'*vY(:) - (wY(:))'*vtY(:))/abs((wtY(:))'*vY(:));
            
            vM = this.DLDm(wM,u,uk,0);
            vY = this.DLDyk(wY,u,uk,0);
            
            vtM = this.DLDm(wtM,u,uk,1);
            vtY = this.DLDyk(wtY,u,uk,1);

            err_DLDm_t = ((wtM(:))'*vM(:) - (wM(:))'*vtM(:))/abs((wtM(:))'*vM(:));
            err_DLDyk_t = ((wtY(:))'*vY(:) - (wY(:))'*vtY(:))/abs((wtY(:))'*vY(:));
            
            results.error_DNDm_t  = err_DNDm_t;
            results.error_DNDy_t  = err_DNDy_t;
            results.error_DLDm_t  = err_DLDm_t;
            results.error_DLDyk_t = err_DLDyk_t;
            
            % =========================================================================

%             report = '';
     
        end
    end

    methods (Access = private)
 
        function set_basis(this)
            if this.rosenbrock == 0 && this.sd == 0
                this.Fourier_basis = 1;
            else
                this.Fourier_basis = 0;
            end
        end

        function Lhatu = Lhat_times(this,u,varargin)

            if this.sd == 0
                Nx = this.grid.x.N;
                Ny = this.grid.y.N;

                uhat  = fft2(reshape(u,Nx,Ny));
                Lhatu = real(ifft2(this.Lhat.*uhat));
                Lhatu = Lhatu(:);
            else
                transpose = 0;
                if nargin > 2
                    transpose = varargin{1};
                end
                if transpose
                    Lhatu = this.Lhat'*u;
                else
                    Lhatu = this.Lhat*u;
                end
            end
        end
        
        function setup_IC(this,varargin)
            
            noise = 1e-1;
            rng(200,'v4');
            this.u0 = noise*randn(this.grid.x.N,this.grid.y.N);
            this.u0 = this.u0(:);

%             xx = linspace(this.grid.x.min,this.grid.x.max,this.grid.x.N+1); xx = xx(1:end-1);
%             yy = linspace(this.grid.y.min,this.grid.y.max,this.grid.y.N+1); yy = yy(1:end-1);
% 
%             [X,Y] = meshgrid(xx,yy);
            if this.Fourier_basis
                this.u0 = this.change_basis(this.u0);
            end
        end
        
        function setup_Lhat(this)
            
            Nx = this.grid.x.N;
            Lx = this.grid.x.max - this.grid.x.min;

            Ny = this.grid.y.N;
            Ly = this.grid.y.max - this.grid.y.min;
            
            if this.sd == 0
                kx = [0:Nx/2-1 -Nx/2:-1]'*2*pi/Lx;   % Wave numbers
                ky = [0:Ny/2-1 -Ny/2:-1]'*2*pi/Ly;

                this.Lhat = zeros(Nx,Ny);
                for nx = 1:Nx
                    for ny = 1:Ny
                        this.Lhat(nx,ny) = -(1-(kx(nx)^2 + ky(ny)^2))^2;
                    end
                end
            else
                xx  = linspace(this.grid.x.min,this.grid.x.max,Nx+1); dx = xx(2)-xx(1);
                yy  = linspace(this.grid.y.min,this.grid.y.max,Ny+1); dy = yy(2)-yy(1);

                D2x = this.derivative_1D(this.sd,Nx,dx);
                D2y = this.derivative_1D(this.sd,Ny,dy);
                
                D2 = kron(D2x,speye(Ny)) + kron(speye(Nx),D2y);

                this.Lhat = -(speye(Nx*Ny) + D2)^2;
            end
        end
        
        function setup_parameters(this,varargin)
            
            pattern = 'cross';
            if nargin > 1 && ischar(varargin{:})
                if strcmpi(varargin{:},'actual')
                    pattern = 'cross';
                elseif strcmpi(varargin{:},'initial guess')
                    pattern = 'smoothed cross';
                end
            end

            Nx = this.grid.x.N;
            Ny = this.grid.y.N;

            r = ones(Nx,Ny);
            g = ones(Nx,Ny);
            
            this.m_bounds.r.min =  0.01; this.m_bounds.r.max = 2.3;
            this.m_bounds.g.min = -1.2;  this.m_bounds.g.max = 1.2;

            if strcmpi(pattern,'uniform')
                r = 0.04*r;
                g = 1*g;
            elseif strcmpi(pattern,'cross')
                r = 2*r;
                r(:,ceil(Ny/4):floor(3*Ny/4)) = 0.04;

                g = -1*g;
                g(ceil(Nx/4):floor(3*Nx/4),:) = 1;
            elseif strcmpi(pattern,'smoothed cross')
                m_grad = this.get_parameter_gradient();
                
                r = 2*r;
                r(:,ceil(Ny/4):floor(3*Ny/4)) = 0.04;

                g = -1*g;
                g(ceil(Nx/4):floor(3*Nx/4),:) = 1;
                
                m = [r(:); g(:)];

                for i = 1:5000
                    m = m + 0.1*m_grad*m;
                end
                
                sets = this.split_parameters(m);
                r = sets{1};
                g = sets{2};
                
            elseif strcmpi(pattern,'cones')
                x = linspace(this.grid.x.min,this.grid.x.max,this.grid.x.N);
                x_mid = 0.5*(this.grid.x.min+this.grid.x.max);
                y = linspace(this.grid.y.min,this.grid.y.max,this.grid.y.N);
                y_mid = 0.5*(this.grid.y.min+this.grid.y.max);
                [X,Y] = meshgrid(x,y);
                
                r = 1/(x_mid*y_mid)*((X-x_mid).^2 + (Y-y_mid).^2);
                r(r < this.m_bounds.r.min) = this.m_bounds.r.min;
                r(r > this.m_bounds.r.max) = this.m_bounds.r.max;
                
                g = 1 - 1/(x_mid*y_mid)*((X-x_mid).^2 + (Y-y_mid).^2);
                g(g < this.m_bounds.g.min) = this.m_bounds.g.min;
                g(g > this.m_bounds.g.max) = this.m_bounds.g.max;
            end
            
            this.m.r = r(:);
            this.m.g = g(:);
            
            this.Nm.r = length(this.m.r);
            this.Nm.g = length(this.m.g);
            
%             figure, ah{1} = axes;
%             fh_pos = get(gcf,'Position');
%                 fh_pos(3) = 400;
%                 fh_pos(4) = 400;
%                 set(gcf,'Position',fh_pos);
%             figure, ah{2} = axes;
%             fh_pos = get(gcf,'Position');
%                 fh_pos(3) = 400;
%                 fh_pos(4) = 400;
%                 set(gcf,'Position',fh_pos);
%             this.display_parameters(ah);
        end
    end
end
