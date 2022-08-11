classdef osp_Hbase_const < solv_prof
    properties (GetAccess=public, SetAccess=private)
        Re
        Sc
        Q
        Gamma
        Dr
    end
    methods
        %% Constructor
        function obj=osp_Hbase_const(mesh_obj,settings,ini_prof)
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
            p.addParameter('Sc',obj.Sc);
            p.addParameter('Q',obj.Q);
            p.addParameter('Gamma',obj.Gamma);
            p.addParameter('Dr',obj.Dr);

            p.addParameter('Nnewton',obj.Nnewton);
            p.addParameter('epsilon',obj.epsilon);
            p.addParameter('min_iter',obj.min_iter);

            p.addParameter('disp_iter',obj.disp_iter);

            p.parse(varargin{:});

            obj.Re=p.Results.Re;
            obj.Sc=p.Results.Sc;
            obj.Q=p.Results.Q;
            obj.Dr=p.Results.Dr;
            obj.Gamma=p.Results.Gamma; % Starting Value of Gamma

            obj.Nnewton=p.Results.Nnewton;
            obj.epsilon=p.Results.epsilon;
            obj.min_iter=p.Results.min_iter;

            obj.disp_iter=p.Results.disp_iter;

        end
        function [res,lhs]=residue_infun(obj,x,lib_value)
            Drr=lib_value.Drr;
            e0rv=lib_value.e0rv;
            % Extracting mesh_obj
            D1M=obj.mesh_obj.D(1);
            D2M=obj.mesh_obj.D(2);
            r01=obj.mesh_obj.pts;
            R=diag(r01);
            N=obj.mesh_obj.N;
            wint=obj.mesh_obj.wint;

            % Extracting x
            G=x(1);
            U0=x(2:N+1);
            H0=x(N+2:2*N+1);

            % Differentiating
            DU0=D1M*U0;
            D2U0=D2M*U0;
            DH0=D1M*H0;
        %     D2N0=D2M*N0;

            % Residue of equation 1
            lhs1=2*pi*(wint*R*U0)-obj.Q;
            if obj.Q==0
                res1=abs(lhs1);
            else
                res1=abs(lhs1)/obj.Q;
            end
            % Residue of equation 2
            Nbar=ones(N,1);
            lhs2=-G*r01+1/obj.Re*DU0+1/obj.Re*r01.*D2U0+obj.Gamma.*r01.*(exp(H0)-Nbar);
            lhs2(end)=DU0(end);
            lhs2(1)=U0(1);
            if any(U0)
                res2=sqrt(lhs2'*diag(wint)*lhs2)/sqrt(U0'*diag(wint)*U0);
            else
                res2=sqrt(lhs2'*diag(wint)*lhs2);
            end
            % Residue of equation 3
            lhs3=1/obj.Sc/obj.Re.*Drr.*DH0-e0rv;
            lhs3(end)=DH0(end); % Unsure why with this commented out, the solver works better.
            res3=sqrt(lhs3'*diag(wint)*lhs3)/sqrt(H0'*diag(wint)*H0);
            % Residue of equation 4
            lhs4=2*(wint*R*exp(H0))-1.;
            res4=abs(lhs4);
            % Combining Residues
            lhs=[lhs1; lhs2;lhs4; lhs3(2:end)];
            res=res1+res2+res3+res4;

        end

        function lib_value= lib_read(obj,y)
            %e0rv is the interpolated <er> vector in the N points, so for r between 1 and 0.
            %Drr is the diffusivity vector between 1 and 0.
            %Der and DDrr are the matrix whose diagonal terms are der/domegapsi (noted Dpsie0r here) and dDrr/domegapsi
            U0=y(2:obj.mesh_obj.N+1);
            Rf=-obj.mesh_obj.D(1)*U0/obj.Dr;
            beta=[-0.103280798025304];
            lib_value.e0rv=beta*Rf;
            lib_value.Drr=[0.0938916345684582];
            lib_value.Der=beta*ones(length(Rf),1);
            lib_value.DDrr=zeros(length(Rf),1);
        end
        function dfdx=dfdx_infun(obj,y,lib_value)
            Drr=lib_value.Drr;
%             e0rv=lib_value.e0rv;
            Der=lib_value.Der;
            DDrr=lib_value.DDrr;
            % Mesh and Differential Matrices
            R=diag(obj.mesh_obj.pts);
            D1M=obj.mesh_obj.D(1);
            D2M=obj.mesh_obj.D(2);
            N=obj.mesh_obj.N;
            wint=obj.mesh_obj.wint;

            % Extracting x
            H0=y(N+2:2*N+1);

            %Differentiating U
            DH0=D1M*H0;

            %matrix formulation to enable the calculus
            Dt22diag=diag(Drr);

            DH0diag=diag(DH0);
            % Setting Jacobian
            dfdx1=[0  2*pi*wint*R zeros(1,N)];

            dfdx2=[-1.*obj.mesh_obj.pts  1/obj.Re.*(D1M+R*D2M)   obj.Gamma.*R*diag(exp(H0))];
            % boundary conditions for u
            dfdx2(1,:)=0;
            dfdx2(N,:)=0;
            dfdx2(1,2)=1.;               %U0(r=1)=0
            dfdx2(N,2:N+1)=D1M(N,:);     %dU0/dr(r=0)=0

            Lnv=-1/obj.Dr*(1/obj.Sc/obj.Re.*(diag(DDrr)*DH0diag*D1M)...
                -diag(Der)*D1M);
            Lnn=1/obj.Sc/obj.Re.*(Dt22diag*D1M);

            dfdx3 =[ zeros(N,1)     Lnv       Lnn] ;
            % boundary conditions for n
    %         dfdx3(1,:)=0;
            dfdx3(N,:)=0;
    %         dfdx3(1,N+2:2*N+1)=1/Sc/Re*Drr(1)*D1M(1,:);
    %         dfdx3(1,N+2)=dfdx3(1,N+2)-e0rv(1);                       %1/Sc/Re*Drr(r=1)*dN0/dr(r=1)-er(r=1)*N0(r=1)=0
            dfdx3(N,N+2:2*N+1)=D1M(N,:);                             %dN0/dr(r=0)=0

            dfdx4=[0  zeros(1,N) 2*wint*R*diag(exp(H0))];

            dfdx=[dfdx1; dfdx2;dfdx4;dfdx3(2:end,:)];
        end
        function obj=get_ini(obj,x0)
            if nargin>1
                obj.x0=x0;
            else
                G=0;
                U0=2/pi*obj.Q*(1-obj.mesh_obj.pts.*obj.mesh_obj.pts);
                H0=zeros(obj.mesh_obj.N,1);

                obj.x0=[G;U0;H0];
            end
        end
    end
end
