classdef LineSearch < handle

    properties (SetAccess = private)

        c_1       = 0.05;
        c_2       = 0.45;
        alpha_max = 1;
        max_iters = 10;

        alphas = [];
        phis   = [];
        dphis  = [];

        alphas_zoom = [];
        phis_zoom   = [];
        dphis_zoom  = [];

        display_data = [];
        
        fcn_value    = [];
        fcn_gradient = [];
        
        % Axes handle
        axes_handle = [];
    end
    
    methods (Access = public)
        function this = LineSearch(varargin)
            if nargin == 1 && isstruct(varargin{1})
                this.c_1       = varargin{1}.c_1;           % Parameters for strong Wolfe conditions
                this.c_2       = varargin{1}.c_2;           %                  
                this.alpha_max = varargin{1}.alpha_max;     % Maximum step length
                this.max_iters = varargin{1}.max_iters;     % Maximum number of line search iterations
            else
                this.set(varargin{:});
            end
            
            this.display_data.alphas     = [];
            this.display_data.phis       = [];
            this.display_data.dphis      = [];
            this.display_data.zoom_start = 0;
        end
        
        % ----------------------- Mutator functions -----------------------
        
        function set(this, varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter, 'value') || strcmpi(parameter, 'fcn_value') || ...
                       strcmpi(parameter, 'f') || strcmpi(parameter, 'function')
                    this.fcn_value = varargin{i+1};
                elseif strcmpi(parameter, 'gradient') || strcmpi(parameter, 'fcn_gradient') || ...
                       strcmpi(parameter, 'grad_f')
                    this.fcn_gradient = varargin{i+1};
                elseif strcmpi(parameter,'c_1') || strcmpi(parameter,'c1')
                    this.c_1 = varargin{i+1};
                elseif strcmpi(parameter,'c_2') || strcmpi(parameter,'c2')
                    this.c_2 = varargin{i+1};
                elseif strcmpi(parameter,'alpha max')
                    this.alpha_max = varargin{i+1};
                elseif strcmpi(parameter,'max iters') || strcmpi(parameter,'max it') || ...
                       strcmpi(parameter,'max iterations')
                    this.max_iters = max(2, varargin{i+1});
                elseif strcmpi(parameter,'axes handle')
                    this.axes_handle = varargin{i+1};
                elseif strcmpi(parameter,'verbose')
                    this.verbose = varargin{i+1};
                end
            end
        end
        
        % ---------------------- Accessor functions -----------------------

        function out = get(this, varargin)
            
            for i = 1:2:nargin-1
                parameter = varargin{i};

                if strcmpi(parameter, 'value') || strcmpi(parameter, 'fcn_value') || ...
                       strcmpi(parameter, 'f') || strcmpi(parameter, 'function')
                    out = this.fcn_value;
                elseif strcmpi(parameter, 'gradient') || strcmpi(parameter, 'fcn_gradient') || ...
                       strcmpi(parameter, 'grad_f')
                    out = this.fcn_gradient;
                elseif strcmpi(parameter,'c_1') || strcmpi(parameter,'c1')
                    out = this.c_1;
                elseif strcmpi(parameter,'c_2') || strcmpi(parameter,'c2')
                    out = this.c_2;
                elseif strcmpi(parameter,'alpha max')
                    out = this.alpha_max;
                elseif strcmpi(parameter,'max iters') || strcmpi(parameter,'max it') || ...
                       strcmpi(parameter,'max iterations')
                    out = this.max_iters;
                elseif strcmpi(parameter,'axes handle')
                    out = this.axes_handle;
                elseif strcmpi(parameter,'verbose')
                    out = this.verbose;
                end
            end
        end
        
        % -----------------------------------------------------------------
        
        function [x,alpha,f,grad_f] = run(this,x,p,alpha,f,grad_f)

            this.reset();
            
            phi  = f;
            dphi = grad_f'*p;
            
            if dphi > 0
                keyboard
            end

            this.alphas(1) = 0;
            this.phis(1)   = phi;
            this.dphis(1)  = dphi;
            this.update_plot(0,phi,dphi);
            
            [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
            this.store(alpha,phi,dphi);
            this.update_plot(alpha,phi,dphi);
            
            % ---------------------------- Scale alpha ----------------------------

            % If alpha is too large, backtrack until an acceptable alpha is found
            while isnan(f) || f > 1.0e10  
                alpha  = 0.1*alpha;
                [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
            end
            
            % ---------------------------------------------------------------------

            for i = 2:this.max_iters

                if this.check_condition(alpha,phi,dphi)
                    break;
                else
                    if dphi > 0
                        [alpha,f,grad_f,phi,dphi] = this.zoom(x,p);
                        break;
                    else
                        if phi > this.phis(i-1)
                            [alpha,f,grad_f,phi,dphi] = this.zoom(x,p);
                            break;
                        else
                            % If current value of phi is less than the previous one, try extrapolating
                            alpha_ext = this.extrapolate(this.alphas(end-1:end), this.phis(end-1:end), this.dphis(end-1:end));
                            if alpha_ext < 0  % If no minimum is found, try twice the value of the current alpha
                                alpha = min(2*alpha,this.alpha_max);
                            else
                                alpha = alpha_ext;
                            end
                        end
                    end
                end

                % Reset alpha if it gets too small
                if alpha < sqrt(eps), alpha = 0.9; end
                
                [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
                this.store(alpha,phi,dphi);
                this.update_plot(alpha,phi,dphi);
            end
            if i == this.max_iters      % If everything else fails, just pick
                                        % the smallest value and be done with it
                [~,ind] = min(this.phis);
                alpha   = this.alphas(ind);
                [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
            end
            x = x + alpha*p;
            this.update_plot();
        end
    end
        
    methods (Access = private)

        function reset(this)
            
            this.alphas = [];
            this.phis   = [];
            this.dphis  = [];

            this.alphas_zoom = [];
            this.phis_zoom   = [];
            this.dphis_zoom  = [];
            
            this.display_data.alphas = [];
            this.display_data.phis   = [];
            this.display_data.dphis  = [];
            this.display_data.zoom_start = 0;
        end

        function alpha = interpolate(this,alphas,phis,dphis)

            hip = pchipd(alphas, phis, dphis);

            for i = 1:size(hip.coefs,1)
                if sum(isnan(hip.coefs(i,:))) || sum(isinf(hip.coefs(i,:))) || sum(~isreal(hip.coefs(i,:)))
                    alpha = -1;
                    return
                end
            end
            
            rts = [];
            for i = 1:size(hip.coefs,1)
                rts_temp = sort(roots(polyder(hip.coefs(i,:)))) + hip.breaks(i);
                if sum(isnan(rts_temp)) || sum(~isreal(rts_temp)) || ...
                        sum(rts_temp < alphas(1)) == length(rts_temp) || ...
                        sum(rts_temp > alphas(2)) == length(rts_temp)
                    rts_temp = [];
                else
                    rts_temp = rts_temp(rts_temp > alphas(1));
                    rts_temp = rts_temp(rts_temp < alphas(2));
                end
                rts = [rts; rts_temp];
            end
            
            if isempty(rts)
                alpha = -1;
            else
                if mod(length(rts),2) == 0, alpha = rts(end-1); else alpha = rts(end); end
                alpha = min(alpha,this.alpha_max);
            end
        end
        
        function alpha = extrapolate(this,alphas,phis,dphis)

            hip = pchipd(alphas, phis, dphis);

            for i = 1:size(hip.coefs,1)
                if sum(isnan(hip.coefs(i,:))) || sum(~isreal(hip.coefs(i,:)))
                    alpha = -1;
                    return
                end
            end
            
            rts = [];
            for i = 1:size(hip.coefs,1)
                rts_temp = sort(roots(polyder(hip.coefs(i,:)))) + hip.breaks(i);
                if sum(isnan(rts_temp)) || sum(~isreal(rts_temp)) || ...
                        sum(rts_temp < alphas(2)) == length(rts_temp)
                    rts_temp = [];
                else
                    rts_temp = rts_temp(rts_temp > alphas(2));
                end
                rts = [rts; rts_temp];
            end
            
            if isempty(rts)
                alpha = -1;
            else
                alpha = min(rts(1),this.alpha_max);
            end
        end
        
        function [phi,dphi,f,grad_f] = evaluate(this,alpha,x,p)

            f      = this.fcn_value(x + alpha*p);    phi  = f;
            grad_f = this.fcn_gradient(x + alpha*p); dphi = grad_f'*p;
        end
        
        function store(this,alpha,phi,dphi,varargin)
            
            if nargin > 4
                i = length(this.alphas_zoom);

                this.alphas_zoom(i+1) = alpha;
                this.phis_zoom(i+1)   = phi;
                this.dphis_zoom(i+1)  = dphi;
            else
                i = length(this.alphas);

                this.alphas(i+1) = alpha;
                this.phis(i+1)   = phi;
                this.dphis(i+1)  = dphi;
            end
        end
        
        function [alpha,f,grad_f,phi,dphi] = zoom(this,x,p)

            if this.alphas(end-1) < this.alphas(end)
                alpha_lo = this.alphas(end-1); alpha_hi = this.alphas(end);
                phi_lo   = this.phis(end-1);   phi_hi   = this.phis(end);
                dphi_lo  = this.dphis(end-1);  dphi_hi  = this.dphis(end);
            else
                alpha_lo = this.alphas(end);   alpha_hi = this.alphas(end-1);
                phi_lo   = this.phis(end);     phi_hi   = this.phis(end-1);
                dphi_lo  = this.dphis(end);    dphi_hi  = this.dphis(end-1);
            end
            
            for i = 1:this.max_iters
                
                if mod(i,2) == 0                            % Give up and try the midpoint
                    alpha = 0.5*(alpha_lo + alpha_hi);
                else
                    alpha = this.interpolate([alpha_lo alpha_hi],[phi_lo phi_hi],[dphi_lo dphi_hi]);
                end
                
                if alpha < sqrt(eps)        % Give up and try the midpoint if alpha gets too small
                    alpha = 0.5*(alpha_lo + alpha_hi);
                end
                
                [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
                this.store(alpha,phi,dphi,'zoom');
                this.update_plot(alpha,phi,dphi,'zoom');
                
                if this.check_condition(alpha,phi,dphi)
                    break;
                else
                    if dphi > 0 || phi > phi_lo
                        alpha_hi = alpha; phi_hi = phi; dphi_hi = dphi;
                    elseif phi > phi_hi || phi_lo > phi 
                        alpha_lo = alpha; phi_lo = phi; dphi_lo = dphi;
                    else
                        alpha_int = this.interpolate([alpha_lo alpha],[phi_lo phi],[dphi_lo dphi]);
                        if alpha_int < 0
                            alpha_lo = alpha; phi_lo = phi; dphi_lo = dphi;
                        else
                            alpha_hi = alpha; phi_hi = phi; dphi_hi = dphi;
                        end
                    end
                end        
            end

            if i == this.max_iters      % If everything else fails, just pick
                                        % the smallest value and be done with it
                [~,ind] = min(this.phis_zoom);
                alpha   = this.alphas_zoom(ind);
                [phi,dphi,f,grad_f] = this.evaluate(alpha,x,p);
                this.store(alpha,phi,dphi,'zoom');
                this.update_plot(alpha,phi,dphi,'zoom');
            end
        end
        
        function update_plot(this,varargin)
            
            if ~isempty(this.axes_handle)
                if nargin > 1
                    alpha = varargin{1};
                    phi   = varargin{2};
                    dphi  = varargin{3};

                    i = length(this.display_data.alphas);

                    this.display_data.alphas(i+1) = alpha;
                    this.display_data.phis(i+1)   = phi;
                    this.display_data.dphis(i+1)  = dphi;

                    if nargin > 4 && this.display_data.zoom_start == 0
                        this.display_data.zoom_start = i+1;
                    end
                    
                    if length(this.display_data.alphas) > 1
                        this.plot_results();
                    end
                else
                    if length(this.display_data.alphas) > 1
                        this.plot_results('final');
                    end
                end
            end
        end

        function plot_results(this,varargin)
            
            ah = this.axes_handle;

            all_alphas = this.display_data.alphas(:);
            all_phis   = this.display_data.phis(:);
            all_dphis  = this.display_data.dphis(:);

            ind      = 1:length(all_alphas);
            ind_zoom = ind;
            i = this.display_data.zoom_start;
            if i > 0
                ind_zoom = ind_zoom(i:end); ind = ind(1:i-1);
            else
                ind_zoom = []; ind = ind;
            end
            ind_zoom = ind_zoom(max(1,end-5):end);
            ind      = ind(max(1,end-5+length(ind_zoom)):end);
            
            min_alpha = min(all_alphas([ind ind_zoom]));
            max_alpha = max(all_alphas([ind ind_zoom]));
            min_phi   = min(all_phis([ind ind_zoom]));
            max_phi   = max(all_phis([ind ind_zoom]));

            alphas_fine = linspace(min_alpha,max_alpha,101);

            delta_alpha = max_alpha - min_alpha;
            delta_phi   = max_phi - min_phi;

            plot(ah,alphas_fine,pchipd(all_alphas([ind ind_zoom]),all_phis([ind ind_zoom]),all_dphis([ind ind_zoom]),alphas_fine(:)),'g');
            hold(ah,'on');
            plot(ah,alphas_fine,this.phis(1) + this.c_1*alphas_fine*this.dphis(1),'k--');
            h1 = plot(ah,all_alphas(ind),all_phis(ind),'o');
            set(h1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4)
            h2 = plot(ah,all_alphas(ind_zoom),all_phis(ind_zoom),'o');
            set(h2,'MarkerEdgeColor','b','MarkerFaceColor','none','MarkerSize',4)
            if nargin > 1       % Special marker for the final iterate
                h3 = plot(ah,all_alphas(end),all_phis(end),'o');
                set(h3,'MarkerEdgeColor','b','MarkerFaceColor','none','MarkerSize',8)
                plot(ah,[all_alphas(end)-0.1*delta_alpha all_alphas(end)+0.1*delta_alpha],...
                    all_phis(end) + [this.c_2*this.dphis(1)*0.1*delta_alpha -this.c_2*this.dphis(1)*0.1*delta_alpha],'r:');
                plot(ah,[all_alphas(end)-0.1*delta_alpha all_alphas(end)+0.1*delta_alpha],...
                    all_phis(end) + [-this.c_2*this.dphis(1)*0.1*delta_alpha this.c_2*this.dphis(1)*0.1*delta_alpha],'r:');
            end

            for i = [ind ind_zoom]

                x = 0.04*delta_alpha;
                y = x*all_dphis(i);

                plot(ah,[all_alphas(i)-x,all_alphas(i)+x],[all_phis(i)-y,all_phis(i)+y],'r');
                
                text(all_alphas(i)+0.01*delta_alpha,all_phis(i)+0.1*delta_phi,num2str(i-1),...
                    'Color', 'k', 'FontSize', 14, 'Fontweight', 'bold','Parent',ah);
            end
            hold(ah, 'off');

            title(ah, {'Cubic spline interpolation using indicated'; 'points (blue disks or circles) and slopes (red line segments)'}, 'FontSize', 10)
            xlabel(ah, '\alpha','FontSize',12)
            ylabel(ah, '\phi (\alpha)','FontSize',12,'Rotation',0)

            drawnow;
        end
          
        function tf = check_condition(this,alpha,phi,dphi)
            tf = ((phi < this.phis(1) + this.c_1*alpha*this.dphis(1)) && ...
                        abs(dphi) <= -this.c_2*this.dphis(1));
        end
    end
end
