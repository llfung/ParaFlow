clear all;
%% ------ Adding paths --------------------------------------------------
addpath('./');
addpath('./core/');
addpath('./core/mesh');
addpath('./core/profile/');
% addpath('./core/stab/');
addpath(genpath('./app/'));
%% ------ Initialization --------------------------------------------------
% Parameters
settings.Re=1;      % Reynolds number based on half width and swimming speed
settings.Gamma=6;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=1;      % Swimming Peclet Number
%settings.Q=0.01;          % Flow rate setting
settings.G0=-0.001;

%% ------ Mesh Initialization ---------------------------------------------
N_mesh=175; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,2,bc_type.none,tran_type.lin);
% mesh=cfd6(N_mesh,2,bc_type.none);

%% ------ Continuation ----------------------------------------------------
settings.disp_iter=true;
settings.epsilon=1e-8;
% Limiter (limiting the value of N(0).
settings.limiter_value=2; % Stop iteration if N(0)>1000
settings.limiter_index=N_mesh+1; % Index for N(0)

U0=zeros(N_mesh,1);%0.00001*cos(pi*mesh.col_pt);
U0=U0-U0(1);
Q=mesh.wint*U0;
H0=zeros(N_mesh,1);

%% Loop Through different starts
Gamma_start=[6 35 85 150];
for i=1:numel(Gamma_start)
    settings.Gamma=Gamma_start(i); 
    cont_obj(i)=NSTransH_VDlib_contRi_G0_sym(mesh,settings,0.05,200,5,[Q;U0;H0]);
    cont_obj(i)=cont_obj(i).run_cont_forward();
end

%% Plotting
figure(4);hold on;
% plot(cont_obj.para_out,cont_obj.x_out(1+floor((N_mesh+1)/2),:));
for i=1:numel(Gamma_start)
    plot(cont_obj(i).para_out,max(cont_obj(i).x_out(2:1+N_mesh,:),[],1));
end


