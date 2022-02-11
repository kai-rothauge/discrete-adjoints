classdef FD_Grid < Grid
    
    properties (SetAccess = protected)
        
        grid_axis = [];         % Grid axis objects
        labels = [];
        
        staggered = 0;
        sd = 2;                 % Spatial discretization order
                               	%   0: pseudospectral
                              	%   2: 2nd
                               	%   4: 4th
                               	%   6: 6th
                                
        n  = [];                % Nodes of grid
        cf = [];                % Cell faces of grid
        cc = [];                % Cell centers of grid
        
        An2ce  = [];
        Acf2cc = [];
        An2cc  = [];
    end
    
    methods (Access = public)
        function this = FD_Grid(varargin)
            
            this.set_default_labels();
            this.initialize_grid_axes();
            if nargin == 1
                this.set('labels',varargin{:});
            else
                this.set(varargin{:});
            end
        end
        
        function flag = set_default_labels(this)
            
            if isempty(this.labels)
                labels = [];
                for i = 1:this.dim
                    labels{i} = ['x' num2str(i)];
                end
                this.set_labels(labels);
            end
            flag = 0;
        end
        
        function flag = map_labels(this)
            
            for i = 1:length(this.grid_axis)
                this.labels{i} = this.grid_axis.get('label');
            end
            flag = 0;
        end
        
        function flag = initialize_grid_axes(this)

            for i = 1:this.dim
                this.grid_axis{i} = Grid_Axis(this.labels(i));
            end
            flag = 0;
        end
        
        % ----------------------- Mutator functions -----------------------
            
        function set(this,varargin)
            
            if nargin > 1
                for i = 1:2:nargin-1
                    parameter = varargin{i};

                    if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                        disp('  ''dim''');
                        disp('  ''labels'' or ''axis labels''');
                        disp('  ''staggered'' or ''staggered grid''');
                    elseif strcmpi(parameter,'dim')
                        this.set_dim(varargin{i+1});
                    elseif strcmpi(parameter,'labels') || strcmpi(parameter,'axis labels')
                        this.set_labels(varargin{i+1});
                    elseif strcmpi(parameter,'staggered') || strcmpi(parameter,'staggered grid')
                        this.set_staggered(varargin{i+1});
                    end
                end
            end
            
            this.rebuild();
        end
        
        function flag = rebuild(this)
            this.build_grid();
            this.build_averaging_operators();
            this.build_derivative_operators();
            flag = 0;
        end
        
        % ----------------------- Accessor functions ----------------------
        
        function out = get(this,varargin)
            
            out = [];
            if nargin > 1
                parameter = varargin{1};

                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('  ''dim''');
                    disp('  ''axis'' or ''grid axis''');
                    disp('  ''staggered'' or ''staggered grid''');
                elseif strcmpi(parameter,'dim')
                    out = this.dim;
                elseif strcmpi(parameter,'axis') || strcmpi(parameter,'grid axis')
                    out = this.dim;
                elseif strcmpi(parameter,'staggered') || strcmpi(parameter,'staggered grid')
                    this.set_staggered(varargin{i+1});
                end
            end
        end
        
    end

    methods (Access = private)
        
        % ----------------------- Mutator functions -----------------------
               
        function set_labels(this,labels)
            this.set_dim(length(labels));
            this.labels = labels;
            if ~isempty(this.grid_axis)
                for i = 1:this.dim
                    this.grid_axis{i}.set('label',labels{i});
                end
            end
        end
        
        function set_staggered(this,staggered)
            this.staggered = (staggered > 0);
        end
        
        function flag = set_dim(this,dim)
            
            flag = 0;
            if dim < 1 || dim > 3, flag = 1;
            else this.dim = dim; end
        end
        
        % -----------------------------------------------------------------

        function build_grid(this)
            
            if this.dim == 1
                
                
            elseif this.dim == 2
                
                xn = this.grid_axis{1}.node;
                xm = this.grid_axis{1}.mid;
                yn = this.grid_axis{2}.node;
                ym = this.grid_axis{2}.mid;
                
                this.n     = zeros(length(xn),length(yn));
                this.cc    = zeros(length(xm),length(ym));
                this.cf{1} = zeros(length(xn),length(ym));
                this.cf{2} = zeros(length(xm)*length(yn));
                
                for i = 1:length(xn)
                    for j = 1:length(yn)
                        this.n(i,j) = Grid_Point(xn(i),yn(j));
                    end
                end
                
                for i = 1:length(xm)
                    for j = 1:length(ym)
                        this.cc(i,j) = Grid_Point(xm(i),ym(j));
                    end
                end
                
                for i = 1:length(xn)
                    for j = 1:length(ym)
                        this.cf{1}(i,j) = Grid_Point(xn(i),ym(j));
                    end
                end
                
                for i = 1:length(xm)
                    for j = 1:length(yn)
                        this.cf{2}(i,j) = Grid_Point(xm(i),yn(j));
                    end
                end
                    
            elseif this.dim == 3
                
            end
            
        end
        
        function build_averaging_operators(this)

            if this.dim == 1
                
                
            elseif this.dim == 2
                
                N1 = this.grid_axis{1}.get('N');
                N2 = this.grid_axis{2}.get('N');

                A1 = 0.5*spdiags(ones(N1,2),[0 1],N1-1,N1);
                A2 = 0.5*spdiags(ones(N2,2),[0 1],N2-1,N2);

                this.An2ce = [kron(A1,speye(N2));
                              kron(speye(N1),A2)];

                this.Acf2cc = [kron(A1, speye(N2-1)) kron(speye(N1-1), A2)];

                this.An2cc = 0.25*spdiags(ones((N1-1)*N2, 4),[0 1 (N2+1) (N2+2)], (N1-1)*N2, N1*N2);
                this.An2cc(N2*(1:N1-1),:) = [];
            elseif this.dim == 3
                
            end
        end

        function build_derivative_operators(this)

            if this.dim == 1
                
                
            elseif this.dim == 2

                D1 = this.grid_axis{1}.get_derivative_operator(this.sd,this.staggered);
                D2 = this.grid_axis{2}.get_derivative_operator(this.sd,this.staggered);

                N1 = this.grid_axis{1}.get('N');
                N2 = this.grid_axis{2}.get('N');
                
                this.Grad = [kron(speye(N2),D1); kron(D2,speye(N1))];
                this.Div  = -this.Grad';
                this.Curl = [kron(D1', speye(N2)) -kron(speye(N1), D2')];
                this.Lap  = -this.Div*this.Grad;
            elseif this.dim == 3
                
            end
        end
        
    end
end
