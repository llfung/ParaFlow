classdef NSTransH_VDlib < solv_prof
    properties (GetAccess=public, SetAccess=private)
        Re
        Pe_s
        Q
        Gamma
    end
    methods
        %% Constructor
        function obj=NSTransH_VDlib(mesh_obj,settings,ini_prof)
            switch nargin
                case 2
                    super_args{1}=mesh_obj;
                    super_args{2}=settings;
                    super_args{3}=2*double(mesh_obj.N)+1; % Size of Prof
                case 3
                    super_args{1}=mesh_obj;
                    super_args{2}=settings;
                    super_args{3}=2*double(mesh_obj.N)+1; % Size of Prof
                    super_args{4}=ini_prof;
                otherwise
                    error('incorrect number of inputs');
            end
            obj=obj@solv_prof(super_args{:});
        end
        %% Inherited Abstract class
        function obj=set_para(obj,varargin)
            narginchk(2,21);
            p=inputParser;
            p.KeepUnmatched=true;
            p.PartialMatching=false;
            p.addParameter('Re',obj.Re);
            p.addParameter('Pe_s',obj.Pe_s);
            p.addParameter('Q',obj.Q);
            p.addParameter('Gamma',obj.Gamma);
            
            p.addParameter('Nnewton',obj.Nnewton);
            p.addParameter('epsilon',obj.epsilon);
            p.addParameter('min_iter',obj.min_iter);

            p.addParameter('disp_iter',obj.disp_iter);


            p.parse(varargin{:});

            obj.Re=p.Results.Re;
            obj.Pe_s=p.Results.Pe_s;
            obj.Q=p.Results.Q;
            obj.Gamma=p.Results.Gamma; 

            obj.Nnewton=p.Results.Nnewton;
            obj.epsilon=p.Results.epsilon;
            obj.min_iter=p.Results.min_iter;

            obj.disp_iter=p.Results.disp_iter;

        end
        function [res,lhs]=residue_infun(obj,x,lib_value)
            D11=lib_value.D11*obj.Pe_s^2;
            V1=lib_value.V1*obj.Pe_s;
            % Extracting mesh_obj
            D1M=obj.mesh_obj.D(1);
            D2M=obj.mesh_obj.D(2);
            N=double(obj.mesh_obj.N);
            wint=obj.mesh_obj.wint;

            % Extracting x
            G=x(1);
            U0=x(2:N+1);
            H0=x(N+2:2*N+1);

            % Differentiating
            DU0=D1M*U0;
            D2U0=D2M*U0;
            DH0=D1M*H0;
%             D2H0=D2M*H0;

            % Residue of equation 1
            lhs1=(wint*U0)-obj.Q*2;
            if obj.Q==0
                res1=abs(lhs1);
            else
                res1=abs(lhs1)/obj.Q/2;
            end
            % Residue of equation 2
            Nbar=ones(N,1);
            lhs2=-G+1/obj.Re*D2U0-obj.Gamma*(exp(H0)-Nbar);
            lhs2(end)=0;
            lhs2(1)=0;
            if any(U0) && sqrt(U0'*diag(wint)*U0)>5e-4
                res2=sqrt(lhs2'*diag(wint)*lhs2)/sqrt(U0'*diag(wint)*U0);
            else
                res2=sqrt(lhs2'*diag(wint)*lhs2);
            end
            lhs2(end)=U0(end);
            lhs2(1)=U0(1);

            % Residue of equation 3
            lhs3=D11.*DH0-V1;
            res3=sqrt(lhs3'*diag(wint)*lhs3)/sqrt(H0'*diag(wint)*H0);
            
            % Residue of equation 4
            lhs4=(wint*exp(H0))-2.;
            res4=abs(lhs4);

            % Combining Residues
            lhs=[lhs1; lhs2; lhs4; lhs3(2:end)];
            res=res1+res2+res3+res4;

        end
        function lib_value = lib_read(obj,y)
            %e0rv is the interpolated <er> vector in the N points, so for r between 1 and 0.
            %Drr is the diffusivity vector between 1 and 0.
            %er and DDrr are the matrix whose diagonal terms are der/domegapsi (noted Dpsie0r here) and dDrr/domegapsi
            persistent loadmat
            if isempty(loadmat)
                loadmat=load('./GTD3D_libp_beta22_fit.mat','e1fitobject','D11fitobject');
            end
            %% Compute S
            U0=y(2:obj.mesh_obj.N+1);
            S=obj.mesh_obj.D(1)*U0;
            %% Interpolate
            lib_value.V1=loadmat.e1fitobject(S);
            lib_value.D11=loadmat.D11fitobject(S);
            lib_value.DV1=differentiate(loadmat.e1fitobject,S);
            lib_value.DD11=differentiate(loadmat.D11fitobject,S); 
        end
        function dfdx=dfdx_infun(obj,y,lib_value)
            V1=lib_value.V1*obj.Pe_s;
            DV1=lib_value.DV1*obj.Pe_s;
            D11=lib_value.D11*obj.Pe_s^2;
            DD11=lib_value.DD11*obj.Pe_s^2;
            % Mesh and Differential Matrices
            D1M=obj.mesh_obj.D(1);
            D2M=obj.mesh_obj.D(2);
            N=obj.mesh_obj.N;
            wint=obj.mesh_obj.wint;

            % Extracting x
            H0=y(N+2:2*N+1);
            DH0=D1M*H0;

            % Setting Jacobian
            dfdx1=[0  wint zeros(1,N)];

            dfdx2=[-ones(obj.mesh_obj.N,1) 1/obj.Re.*(D2M)   -obj.Gamma*diag(exp(H0))];
            % boundary conditions for u
            dfdx2(1,:)=0;
            dfdx2(N,:)=0;
            dfdx2(1,2) = 1.;                 % U0(x = 1)=0
            dfdx2(N,N+1) = 1;     % U0(x =-1)=0

            Lnv=diag(DH0.*DD11)*D1M...
                -diag(DV1)*D1M;
            Lnn=diag(D11)*D1M;

            dfdx3 =[ zeros(N,1)     Lnv       Lnn] ;

            dfdx4=[0  zeros(1,N) wint*diag(exp(H0))];

            dfdx=[dfdx1; dfdx2;dfdx4;dfdx3(2:end,:)];
        end
        function obj=get_ini(obj,x0)
            if nargin>1
                obj.x0=x0;
            else
                G=0;
                U0=obj.Q/2*3*((1-obj.mesh_obj.pts.^2));%+sin(pi*obj.mesh_obj.pts));
                H0=zeros(obj.mesh_obj.N,1);

                obj.x0=[G;U0;H0];
            end
        end
    end
end
