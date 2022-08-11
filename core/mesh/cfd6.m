classdef cfd6 < mesh_class
%     Properties inherented from mesh_class
%         pts %Mesh Points
%         DM % Differential Matrices (3rd dimenion is the order of differentiation)
%         N uint16 %Size of Matrices and Collocation Points Vector
%         M uint8 %Max Order of Chebyshev Differential Matrices
%         bc_opt bc_type %Boundary Type
%         wint %Integration weight vector
    methods
        %% Constructor
        function obj=cfd6(siz,M,bc_opt)
            narginchk(1,3);
            if nargin < 2
                M=2;
            elseif ~(M>=0 && M <=2) 
                error('M value incorrect. Must be integer between 0 to 2.');
            end
            obj=obj@mesh_class(siz,M);

            %% Forming the differential matrix
            if obj.M >0
                obj.DM(:,:,1)=spdiags(ones(obj.N,1)*[-1/60 3/20 -3/4 0 3/4 -3/20 1/60],[-3:1:3],obj.N,obj.N);
            end
            if obj.M >1
                obj.DM(:,:,2)=spdiags(ones(obj.N,1)*[1/90 -3/20 3/2 -49/18 3/2 -3/20 1/90],[-3:1:3],obj.N,obj.N);
            end
            %% Boundary Conditions
            if nargin < 3
                obj.bc_opt=bc_type.none;
            end
            switch bc_opt
                case bc_type.none
                    d=2/double(obj.N-1);
                    obj.pts = (-1:d:1)';
                    if obj.M >0
                        % TODO: Higher Degree FD
                        obj.DM(1:3,:,1)=0;
                        obj.DM(1,1:4,1)=[-11/6 3 -3/2 1/3];
                        obj.DM(2,1:3,1)=[-1/2 0 1/2];
                        obj.DM(3,1:5,1)=[1/12 -2/3 0 2/3 -1/12];

                        obj.DM(end-2:end,:,1)=0;
                        obj.DM(end,end-3:end,1)=[-1/3 3/2 -3 11/6];
                        obj.DM(end-1,end-2:end,1)=[-1/2 0 1/2];
                        obj.DM(end-2,end-4:end,1)=[1/12 -2/3 0 2/3 -1/12];
                    end
                    if obj.M >1
                        % TODO: Higher Degree FD
                        obj.DM(1:3,:,2)=0;
                        obj.DM(1,1:4,2)=[2 -5 4 -1];
                        obj.DM(2,1:3,2)=[1 -2 1];
                        obj.DM(3,1:5,2)=[-1/12 4/3 -5/2 4/3 -1/12];

                        obj.DM(end-2:end,:,2)=0;
                        obj.DM(end,end-3:end,2)=[-1 4 -5 2];
                        obj.DM(end-1,end-2:end,2)=[1 -2 1];
                        obj.DM(end-2,end-4:end,2)=[-1/12 4/3 -5/2 4/3 -1/12];
                    end
                    obj.wint=ones(1,obj.N)*d;
                    obj.wint(1)=obj.wint(1)/2;
                    obj.wint(end)=obj.wint(end)/2;
                case bc_type.periodic
                    d=2/double(obj.N);
                    obj.pts = (-1:d:1-d)';
                    if obj.M >0
                        obj.DM(:,:,1)=spdiags(ones(obj.N,1)*[-1/60 3/20 -3/4 3/4 -3/20 1/60],[double(obj.N)-3:1:double(obj.N)-1 -double(obj.N)+1:1:-double(obj.N)+3],obj.DM(:,:,1));
                    end
                    if obj.M >1
                        obj.DM(:,:,2)=spdiags(ones(obj.N,1)*[1/90 -3/20 3/2 3/2 -3/20 1/90],[double(obj.N)-3:1:double(obj.N)-1 -double(obj.N)+1:1:-double(obj.N)+3],obj.DM(:,:,2));
                    end
                    obj.wint=ones(1,obj.N)*d;
                otherwise
                    error('Other boundary condition types have yet to be implemented.');
            end
            obj.DM(:,:,1)=obj.DM(:,:,1)/d;
            obj.DM(:,:,2)=obj.DM(:,:,2)/d^2;
        end
        %% Inherented abstract class
        function integ=int(obj,v)
            integ=sum(v)/double(obj.N);
        end
    end
end

