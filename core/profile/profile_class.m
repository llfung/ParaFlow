classdef profile_class
    properties
        mesh_obj
        prof %Profile Vector 
        N %Size of profile vector
    end
    methods (Abstract)
        set_prof(obj,settings)
    end
    methods
        function obj=profile_class(mesh_obj,N)
            obj.mesh_obj=mesh_obj;
            if nargin >2
                obj.N=N;
            else
                obj.N=mesh_obj.N;
            end
            obj.prof=zeros(obj.N,1);
        end
    end
end
