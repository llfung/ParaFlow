classdef solv_prof < profile_class
    properties (GetAccess=public, SetAccess=protected)   
        res
        n
        x0
        
        Nnewton(1,1) double {mustBePositive,mustBeInteger} = 1200
        epsilon(1,1) double {mustBePositive} = 2.5e-10
        min_iter(1,1) double {mustBePositive,mustBeInteger} = 50
        
        disp_iter(1,1) logical = false
    end
    methods (Abstract)
        [res,lhs]=residue_infun(obj,y,lib_value)
        dfdx=dfdx_infun(obj,y,lib_value)
        lib_value=lib_read(obj,y)
        obj=set_para(obj,varargin)
        obj=get_ini(obj,x0)
    end
    methods
        function obj=solv_prof(mesh_obj,settings,N_prof,ini_prof)
            obj=obj@profile_class(mesh_obj,N_prof);
            obj=obj.set_para(settings);
            if nargin>3
                if ischar(ini_prof) || isstring(ini_prof)
                    load([ini_prof '.mat'],'x0');
                    obj=obj.get_ini(x0);
                else
                    obj=obj.get_ini(ini_prof);
                end
            else
                obj=obj.get_ini();
            end
        end
        function obj=set_prof(obj,settings)
            if nargin >1
                obj=obj.set_para(settings);
                if isfield(settings,'x0')
                    obj=obj.get_ini(settings.x0);
                elseif isempty(obj.prof)
                    obj=obj.get_ini();
                else
                    obj=obj.get_ini(obj.prof);
                end
            end
            obj=solve(obj);
        end
        function obj=solve(obj,x0)
            %% Initialisation
            if nargin > 1
                obj.x0=x0;
            elseif isempty(obj.x0)
                obj=obj.get_ini();
            end
            x=obj.x0;
            
            res=Inf;
            i=1;
            %% Newton-Raphson Method loop
            while (i<=obj.Nnewton && res>obj.epsilon)  || i< obj.min_iter
                try
                    lib_value=obj.lib_read(x);
                catch
                    break;
                end

                [res,lhs]=obj.residue_infun(x,lib_value);

                dfdx=obj.dfdx_infun(x,lib_value);

                ls=dfdx\lhs;

                x=x-ls;
                i=i+1;
            end
            obj.res=res;
            obj.n=i-1;
            obj.prof=x;
            if obj.disp_iter
                disp(['Iteration finished in ' num2str(obj.n) ' iterations with residue value: ' num2str(obj.res)]);
            end
        end

    end
end
  