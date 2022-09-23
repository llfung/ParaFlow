classdef NSTransH_VDconst_radial < solv_prof
    properties (GetAccess=public, SetAccess=private)
        Re
        Pe_s
        Q
        Gamma
    end
    methods
        %% Constructor
        function obj=NSTransH_VDconst_radial(mesh_obj,settings,ini_prof)
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
            rinv=1./obj.mesh_obj.col_pt;rinv(1)=0;
            wint=obj.mesh_obj.wint.*obj.mesh_obj.col_pt';

            % Extracting x
            G=x(1);
            U0=x(2:N+1);
            H0=x(N+2:2*N+1);

            % Differentiating
            DU0=D1M*U0;
            D2U0=D2M*U0;
            DH0=D1M*H0;
            % D2H0=D2M*H0;

            % Residue of equation 1 (Flow Rate)
            lhs1=2*(wint*U0)-obj.Q;
            if obj.Q==0
                res1=abs(lhs1);
            else
                res1=abs(lhs1)/obj.Q;
            end

            % Residue of equation 2 (Navier-Stokes)
            Nbar=ones(N,1);
            lhs2=-G+1/obj.Re*(rinv.*DU0+D2U0)-obj.Gamma.*(exp(H0)-Nbar);
            lhs2(end)=0;
            lhs2(1)=0;
            if any(U0) && sqrt(U0'*diag(wint)*U0)>5e-4
                res2=sqrt(lhs2'*diag(wint)*lhs2)/sqrt(U0'*diag(wint)*U0);
            else
                res2=sqrt(lhs2'*diag(wint)*lhs2);
            end
            lhs2(end)=DU0(end);
            lhs2(1)=U0(1);

            % Residue of equation 3 (Transport, integrated)
            lhs3=D11.*DH0-V1;
            lhs3(end)=0;
            res3=sqrt(lhs3'*diag(wint)*lhs3)/sqrt(H0'*diag(wint)*H0);
            lhs3(end)=DH0(end);
            
            % Residue of equation 4 (Conservation of N)
            lhs4=2*(wint*exp(H0))-1.;
            res4=abs(lhs4);

            % Combining Residues
            lhs=[lhs1; lhs2; lhs4; lhs3(2:end)];
            res=res1+res2+res3+res4;

        end
        function lib_value = lib_read(obj,y)
            %V1 is drift in r direction.
            %D11 is the diffusivity in r direction.
            %DV1 and DD11 are dV1/dS and dD11/dS
            %% Compute S
            U0=y(2:obj.mesh_obj.N+1);
            S=obj.mesh_obj.D(1)*U0;
            %% Interpolate
            lib_value.V1=-S;
            lib_value.D11=ones(obj.mesh_obj.N,1);
            lib_value.DV1=-ones(obj.mesh_obj.N,1);
            lib_value.DD11=zeros(obj.mesh_obj.N,1);
        end
        function dfdx=dfdx_infun(obj,y,lib_value)
            %V1=lib_value.V1*obj.Pe_s;
            DV1=lib_value.DV1*obj.Pe_s;
            D11=lib_value.D11*obj.Pe_s^2;
            DD11=lib_value.DD11*obj.Pe_s^2;
            % Mesh and Differential Matrices
            D1M=obj.mesh_obj.D(1);
            D2M=obj.mesh_obj.D(2);
            N=obj.mesh_obj.N;
            rinv=1./obj.mesh_obj.col_pt;rinv(1)=0;
            wint=obj.mesh_obj.wint.*obj.mesh_obj.col_pt';

            % Extracting x
            H0=y(N+2:2*N+1);
            DH0=D1M*H0;

            % Preparing Jacobian
            % (Flow Rate)
            dfdx1=[0  2*wint zeros(1,N)];

            % (Navier-Stokes)
            dfdx2=[-ones(obj.mesh_obj.N,1)  1/obj.Re*(diag(rinv)*D1M+D2M)  -obj.Gamma*diag(exp(H0))];
            % boundary conditions for u
            dfdx2(1,:) = 0.;
            dfdx2(N,:) = 0.;
            dfdx2(1,2) = 1.;                 %  U0(x = 1)=0
            dfdx2(N,2:N+1) = D1M(end,:);     % DU0(x = 0)=0
            
            % (Transport)
            Lnv=diag(DH0.*DD11)*D1M...
                -diag(DV1)*D1M;
            Lnn=diag(D11)*D1M;

            dfdx3 =[ zeros(N,1)     Lnv       Lnn] ;
            % boundary conditions for N
            dfdx3(N,:) = 0;
            dfdx3(N,N+2:2*N+1) = D1M(N,:);     % DH0(x = 0)=0
            
            % (Conservation of N)
            dfdx4=[0  zeros(1,N) 2*wint*(diag(exp(H0)))];

            % Combining Jacobians
            dfdx=[dfdx1; dfdx2; dfdx4; dfdx3(2:end,:)];
        end
        function obj=get_ini(obj,x0)
            if nargin>1
                obj.x0=x0;
            else
                G=0;
                U0=obj.Q*4*(1-obj.mesh_obj.pts.^2);
                H0=zeros(obj.mesh_obj.N,1);

                obj.x0=[G;U0;H0];
            end
        end
    end
end
