classdef (Abstract) ode_continuation
    %ode_continuation Object Class for continuation algorithm of an 1D ODE
    %   TODO: 
    %     - Description of this class
    %     - made as inherent class from prof_class
    %     - Support of getting profile solution similar to solv_prof for
    %     stability analysis downstream

    
    properties (GetAccess=public, SetAccess=protected)
        %% Inputs
        mesh_obj %chebyshev
        x0 %Initial State (Vector)
        
        %% Options (Madatary ones)
        % Iterations
        num_iter  (1,1) double {mustBePositive,mustBeInteger} = 1  % Number of iteration saved 
        save_iter (1,1) double {mustBePositive,mustBeInteger} = 1  % Frequency of saving
        total_iter(1,1) double {mustBePositive,mustBeInteger} = 1  % Number of iteration = num_iter * save_iter
        % Continuation step size
        h 
        % Stopping Criteria
        Nnewton(1,1) double {mustBePositive,mustBeInteger} = 3000
        epsilon(1,1) double {mustBePositive} = 1e-9
        min_iter(1,1) double {mustBePositive,mustBeInteger} =3
        % Limiter stopping criteria
        limiter_index
        limiter_value
        %% Parameters
        iter_para string
        
        %% Output
        x_out
        para_out
        res_out
        n_out

        %% Display
        disp_iter(1,1) logical = false
    end
    %% Abstract Classes
    methods (Abstract)
        [res,lhs]=residue_infun(obj,y,yp,lib_value,z)
        dfdx=dfdx_infun(obj,y,lib_value,z)
        lib_value=lib_read(obj,y)
        obj=set_para(obj,varargin)
        x0=get_initial_solution(obj)
    end
    methods
        %% Constructor
        function obj = ode_continuation(mesh_obj,settings,iter_para,h,num_iter,save_iter,x0)
            % ode_continuation Construct an instance of this class
            %   Detailed explanation goes here
            
            %% Initialization
            obj.mesh_obj=mesh_obj;
            obj.h=h;
            obj.num_iter=num_iter;
            obj.iter_para=iter_para;
            if nargin>5
                obj.save_iter=save_iter;
            else
                obj.save_iter=1;
            end
            obj.total_iter=obj.save_iter*obj.num_iter;
            obj=obj.set_para(settings);

            % Put in initial condition
            if nargin>6 && ~isempty(x0)
                obj.x0=x0;
            else
                obj.x0=obj.get_initial_solution();
                if any(isnan(obj.x0))
                    error('ode_continuation: NaN value detected at initial condition');
                end
            end
        end
        
        function obj=run_cont_forward(obj)
            [obj.x_out,obj.para_out,obj.res_out,obj.n_out] = run_cont(obj,obj.h,obj.x0);
        end
        function obj=run_cont_backward(obj)
            [obj.x_out,obj.para_out,obj.res_out,obj.n_out] = run_cont(obj,-obj.h,obj.x0);
        end
        function obj=run_cont_both(obj)
            obj.num_iter=floor(obj.num_iter/2);
            obj.total_iter=obj.save_iter*obj.num_iter;
            [x_out_for,para_out_for,res_out_for,n_out_for] = run_cont(obj,obj.h,obj.x0);
            [x_out_bak,para_out_bak,res_out_bak,n_out_bak] = run_cont(obj,-obj.h,obj.x0);
            obj.x_out   =[    x_out_bak(:,end:-1:1)    x_out_for];
            obj.para_out=[ para_out_bak(1,end:-1:1) para_out_for];
            obj.res_out =[  res_out_bak(1,end:-1:1)  res_out_for];            
            obj.n_out   =[    n_out_bak(1,end:-1:1)    n_out_for];
            obj.num_iter=obj.num_iter*2;
            obj.total_iter=obj.total_iter*2;
        end
        
        function [x_out,para_out,res_out,n_out] = run_cont(obj,h,x0)
            %run_cont Run the Continuation
            %   TODO: Detailed explanation goes here
            
            %% Initialisation
            N_vec=size(x0,1);
            x_out   =NaN(N_vec,obj.num_iter);
            para_out=NaN(1,obj.num_iter);
            n_out   =NaN(1,obj.num_iter);
            res_out =NaN(1,obj.num_iter); 

            iter_para_ini_value=obj.(obj.iter_para);
            
            % Initial extrapolation vector
            z=[zeros(1,N_vec) 1];
            % Initial solution with iterated parameter at the end
            y=[x0;iter_para_ini_value];
            
            breakflag=false;
            
            %warning supression
            warning('off','all');

            %% Iteration for continuation
            for j=1:obj.total_iter
                % Extrapolation with extrapolation vector
                yp=y+h*z'; 
                y=yp;

                res=Inf;
                i=1;
                %% Newton-Raphson Method loop
                while (i<=obj.Nnewton && res>obj.epsilon) || (i<obj.min_iter)
                    try 
                        lib_value=obj.lib_read(y);
                    catch
                         breakflag=true;
                         break;
                    end
                    
                    [res,lhs]=obj.residue_infun(y,yp,lib_value,z);

                    dfdx=obj.dfdx_infun(y,lib_value,z);

                    ls=dfdx\lhs;

                    y=y-ls;
                    i=i+1;
                end

                %% Saving and post-processing after loop
                if mod(j,obj.save_iter)==0
                    if obj.disp_iter
                      disp(strcat(num2str(j),"th profile found in ",num2str(i-1)," iterations with residue value: ",num2str(res),...
                          " with  ",obj.iter_para,"=",num2str(y(end,1))));
                    end
                    x_out   (:,j/obj.save_iter)=y(1:N_vec,1);
                    para_out(1,j/obj.save_iter)=y(N_vec+1,1);
                    n_out   (1,j/obj.save_iter)=i-1;
                    res_out (1,j/obj.save_iter)=res;
                end
                if  res>obj.epsilon || any(isnan(y)) || breakflag
                    break;
                end
                if ~isempty(obj.limiter_index)
                    if abs(y(obj.limiter_index))>obj.limiter_value
                        break;
                    end
                end

                %% Next extrapolation vector
                z=(dfdx\[zeros(N_vec,1);1])';

            end
            %warning back on
            warning('on','all');
        end
    end
end

