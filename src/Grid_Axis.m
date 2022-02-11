classdef Grid_Axis < handle
    
    properties (SetAccess = private)
        
        label = 'x';                % Axis label
        min   = 0;                  % Grid axis min
        max   = 1;                  % Grid axis max
        N     = 32;                 % Number of grid points along axis
        node  = [];                 % Grid points (nodes)
        mid   = [];                 % Midpoints between nodes
        h     = [];                 % Intervals between nodes
        BCs   = ['P' 'P'];          % Boundary conditions
                                    %   D: Dirichlet, N: Neumann
                                    %   P: Periodic,  A: Absorbing
    end
    
    methods (Access = public)
        function this = Grid_Axis(varargin)

            if nargin == 1
                this.set('label',varargin{1});
            else
                this.set(varargin{:});
            end
        end
        
        % ----------------------- Mutator functions -----------------------
          
        function set(this,varargin)
            
            if nargin > 1
                for i = 1:2:nargin-1
                    parameter = varargin{i};

                    if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                        disp('  ''min''');
                        disp('  ''max''');
                        disp('  ''N''');
                        disp('  ''BCs'' or ''boundary conditions''');
                        disp('  ''label'' or ''axis label''');
                    elseif strcmpi(parameter,'min')
                        this.set_min(varargin{i+1});
                    elseif strcmpi(parameter,'max')
                        this.set_max(varargin{i+1});
                    elseif strcmpi(parameter,'N')
                        this.set_N(varargin{i+1});
                    elseif strcmpi(parameter,'BCs') || strcmpi(parameter,'boundary conditions')
                        this.set_BCs(varargin{i+1});
                    elseif strcmpi(parameter,'label') || strcmpi(parameter,'axis label')
                        this.set_label(varargin{i+1});
                    end
                end
            end
            
            this.set_grid_points();
        end
        
        function out = get(this,parameter)
            
            if nargin > 1
                if strcmpi(parameter,'what') || strcmpi(parameter,'what?')
                    disp('  ''min''');
                    disp('  ''max''');
                    disp('  ''N''');
                    disp('  ''BCs'' or ''boundary conditions''');
                    disp('  ''label'' or ''axis label''');
                elseif strcmpi(parameter,'min')
                    out = this.min;
                elseif strcmpi(parameter,'max')
                    out = this.max;
                elseif strcmpi(parameter,'N')
                    out = this.N;
                elseif strcmpi(parameter,'BCs') || strcmpi(parameter,'boundary conditions')
                    out = this.BCs;
                elseif strcmpi(parameter,'label') || strcmpi(parameter,'axis label')
                    out = this.label;
                elseif strcmpi(parameter,'node') || strcmpi(parameter,'nodal')
                    out = this.node;
                elseif strcmpi(parameter,'mid') || strcmpi(parameter,'midpoints')
                    out = this.mid;
                elseif strcmpi(parameter,'dx') || strcmpi(parameter,'h')
                    out = this.h;
                end
            end
        end
        
        % ----------------------- Accessor functions ----------------------

        function D = get_derivative_operator(this,sd,varargin)
            
            if nargin > 2 && varargin{1} > 0
                D = this.get_derivative_operator_staggered(sd);
            else
                D = this.get_derivative_operator_nonstaggered(sd);
            end
        end
    end

    methods (Access = private)
        
        % ----------------------- Mutator functions -----------------------
        
        function flag = set_grid_points(this)
            
            this.node = linspace(this.min,this.max,this.N+1)';
            this.mid  = 0.5*(this.node(2:end) + this.node(1:end-1));
            this.h    = this.node(2:end) - this.node(1:end-1);
            
            if this.BCs(1) == 'P' || this.BCs(2) == 'P'
                this.node = this.node(1:end-1);
                this.BCs(1) = 'P';
                this.BCs(2) = 'P';
            elseif this.BCs(1) == 'D'
                
            elseif this.BCs(2) == 'D'
                
            
            end
            flag = 0;
            
        end
        
        function flag = set_min(this,min)
            
            this.min = min;
            flag = 0;
        end
        
        function flag = set_max(this,max)
            
            this.max = max;
            flag = 0;
        end
        
        function flag = set_N(this,N)
            
            if N > 1
                this.N = N;
                flag = 0;
            else
                flag = 1;
            end
        end
        
        function flag = set_label(this,label)
            
            this.label = label;
            flag = 0;
        end
        
        % -------------------------- Derivatives --------------------------
        
        function D = get_derivative_operator_staggered(this,sd)

            x = this.mid;
            
            if (isempty(x) || (length(x) == 1) || (min(x) < sqrt(eps)))
                D = [];
                return;
            end
            
            N   = this.N;
            h   = this.h;
            BCs = this.BCs;

            if BCs(1) == 'P'
                e = ones(N,1);
                
                D = sparse(N,N);

                A = ones(sd,sd);
                b = zeros(sd,1); b(2) = 1;
                for n = 1:N
                    ind_1 = mod(n-sd/2:n-2,N)+1;
                    ind_2 = mod(n:n+sd/2-2,N)+1;

                    temp = zeros(1, sd);
                    for i = 1:sd/2-1
                        temp(i)      = sum(h(ind_1));
                        temp(sd-i+1) = sum(h(ind_2));

                        ind_1 = ind_1(2:end);
                        ind_2 = ind_2(1:end-1);
                    end
                    temp          =  temp+h(n)/2;
                    temp(1:sd/2) = -temp(1:sd/2);

                    for k = 2:sd
                        A(k,:) = temp.^(k-1);
                    end

                    ind = mod(n-sd/2:n+sd/2-1, N)+1;
                    D(n, ind) = (A\b)';
                end
                
            elseif BCs(dim) == 'D'
                N = length(h);
                D = sparse(N,N+1);

                A = ones(sd+1, sd+1);
                b = zeros(sd+1, 1); b(2) = 1;
                for n = 1:(sd/2-1)

                    temp = cumsum(h(n:-1:1));
                    temp = [-temp(end:-1:1); 0; cumsum(h(n+1:sd))];

                    for k = 1:sd
                        A(k+1,:) = (1/factorial(k))*(temp+h(n)/2)'.^k;
                    end

                    D(n, 1:sd+1) = (A\b)';

                end
                A = ones(sd, sd);
                b = zeros(sd, 1); b(2) = 1;
                for n = sd/2:N-sd/2+1

                    temp = zeros(1, sd);
                    for i = 1:sd/2-1
                        temp(i)           = sum(h(n-sd/2+i:n-1));
%                         temp(sd-i+1) = sum(h(n+1:n+sd/2-i));
                    end
                    temp         =  temp+h(n)/2;
                    temp(1:sd/2) = -temp(1:sd/2);

                    for k = 2:sd
                        A(k,:) = temp.^(k-1);
                    end

                    D(n, n-sd/2+1:n+sd/2) = (A\b)';

                end
                A = ones(sd+1, sd+1);
                b = zeros(sd+1, 1); b(2) = 1;
                for n = N:-1:N+2-sd/2

                    temp = cumsum(h(n:-1:N+1-sd));
                    temp = [-temp(end:-1:1); 0; cumsum(h(n+1:end))];

                    for k = 1:sd
                        A(k+1,:) = (1/factorial(k))*(temp+h(n)/2).^k;
                    end

                    D(n, end-sd:end) = (A\b)';
                end
            end
        end
        
        function D = get_derivative_operator_nonstaggered(this,sd)

            x = this.node;
            
            if isempty(x) || length(x) == 1
                D = [];
                return;
            end
            
            h   = this.h;
            N   = this.N;
            BCs = this.BCs;
            
            if BCs(1) == 'P'
                D    = sparse(N,N);
                A    = ones(sd,sd);
                b    = zeros(sd,1); b(2) = 1;
                temp = zeros(1,sd);
                
                ind = [N-sd/2:N 1:sd/2]';

                for n = 1:N
                    ind = [ind; mod(ind(end)+1,N)]; ind = ind(2:end);
                    if ind(end) == 0, ind(end) = N; end
                    
                    hh = h(ind);
                    
                    for i = 1:sd
                        if i <= sd/2
                            temp(i) = -sum(hh(i:sd/2));
                        else
                            temp(i) = sum(hh(sd/2+1:i));
                        end
                    end

                    for k = 2:sd
                        A(k,:) = temp.^(k-1);
                    end

                    c = (A\b)';
                    D(n,ind) = [c(1:sd/2)'; 0; c(sd/2+1:end)'];
                end
            end
            
%             if sd == 2       % 2nd order by default
%                 D = spdiags([-1*e 0*e 1*e],-1:1,N+1,N+1);
%                 D = D/2;
%             elseif sd == 4
%                 D = spdiags([1*e -8*e 0*e 8*e -1*e],-2:2,N+1,N+1);
%                 D = D/12;
%             elseif sd == 6
%                 D = spdiags([-1*e 9*e -45*e 0*e 45*e -9*e 1*e],-3:3,N+1,N+1);
%                 D = D/60;
%             elseif sd == 8
%                 D = spdiags([3*e -32*e 168*e -672*e 0*e 672*e -168*e 32*e -3*e],-4:4,N+1,N+1);
%                 D = D/840;
%             end
% 
%             if BCs(1) == 'P'
%                 e = ones(N,1);
%                 
%                 if sd == 2       % 2nd order by default
%                     D = spdiags([e -2*e e],-1:1,N,N);
%                     D(1,end) = 1; D(end,1) = 1;
%                 elseif sd == 4
%                     D = spdiags([-e 16*e -30*e 16*e -1*e],-2:2,N,N);
%                     D(1,end-1:end) = [-1 16]; D(2,end)   = -1;
%                     D(end,1:2)     = [16 -1]; D(end-1,1) = -1;
%                     D = D/12;
%                 elseif sd == 6
%                     D = spdiags([2*e -27*e 270*e -490*e 270*e -27*e 2*e],-3:3,N,N);
%                     D(1,end-2:end) = [2 -27 270]; D(2,end-1:end) = [2 -27]; D(3,end) = 2;
%                     D(end,1:3)     = [270 -27 2]; D(end-1,1:2)   = [-27 2]; D(end-2,1) = 2;
%                     D = D/180;
%                 end
% %                 D = D/dx^2;
%             elseif BCs(1) == 'D'
%                 
%                 h = this.x.node(2:end) - this.x.node(1:end-1);
%                 D = sparse(N,N+1);
%                 
%                 A = ones(sd+1, sd+1);
%                 b = zeros(sd+1, 1); b(2) = 1;
%                 
%                 keyboard
%                 for n = 1:(sd/2-1)
% 
%                     temp = cumsum(h(n:-1:1));
%                     temp = [-temp(end:-1:1); 0; cumsum(h(n+1:sd))];
% 
%                     for k = 1:sd
%                         A(k+1,:) = (1/factorial(k))*(temp+h(n)/2)'.^k;
%                     end
% 
%                     D(n, 1:sd+1) = (A\b)';
% 
%                 end
%                 A = ones(sd, sd);
%                 b = zeros(sd, 1); b(2) = 1;
%                 for n = sd/2:N-sd/2+1
% 
%                     temp = zeros(1, sd);
%                     for i = 1:sd/2-1
%                         temp(i)           = sum(h(n-sd/2+i:n-1));
% %                         temp(sd-i+1) = sum(h(n+1:n+sd/2-i));
%                     end
%                     temp         =  temp+h(n)/2;
%                     temp(1:sd/2) = -temp(1:sd/2);
% 
%                     for k = 2:sd
%                         A(k,:) = temp.^(k-1);
%                     end
% 
%                     D(n, n-sd/2+1:n+sd/2) = (A\b)';
% 
%                 end
%                 A = ones(sd+1, sd+1);
%                 b = zeros(sd+1, 1); b(2) = 1;
%                 for n = N:-1:N+2-sd/2
% 
%                     temp = cumsum(h(n:-1:N+1-sd));
%                     temp = [-temp(end:-1:1); 0; cumsum(h(n+1:end))];
% 
%                     for k = 1:sd
%                         A(k+1,:) = (1/factorial(k))*(temp+h(n)/2).^k;
%                     end
% 
%                     D(n, end-sd:end) = (A\b)';
%                 end
        end
    end
end
