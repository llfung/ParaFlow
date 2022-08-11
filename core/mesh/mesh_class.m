classdef mesh_class
    properties
        pts % Mesh Points
        DM  % Differential Matrices (3rd dimenion is the order of differentiation)
        N uint16 %Size of Matrices and Collocation Points Vector
        M uint8 %Max Order of Chebyshev Differential Matrices
        bc_opt bc_type %Boundary Type
        wint %Integration weight vector
    end
    methods (Abstract)
        integ=int(obj,v)
    end
    methods
        %% Constructor
        function obj=mesh_class(N,M)
            obj.N=uint16(N);
            obj.M=uint8(M);
        end
        %% Mesh Class functions
        function D_out=D(obj,M)
            if M>obj.M
                error('Exceed Order of Differentiation available');
            end
            D_out=obj.DM(:,:,M);
        end
        function I_out=I(obj)
            I_out=eye(obj.N);
        end
        function O_Out=O(obj)
            O_Out=zeros(obj.N,obj.N);
        end
        %% Legacy Support
        function pts_legacy=col_pt(obj)
            pts_legacy=obj.pts;
        end
    end
end

