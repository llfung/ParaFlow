clear all;
%% ------ Adding paths --------------------------------------------------
addpath('./');
addpath('./core/');
addpath('./core/mesh');
addpath('./core/profile/');
addpath('./core/stab/');
%% ------ Initialization --------------------------------------------------
% Parameters
settings.Re=0.126;       % Reynolds number based on half width and swimming speed
settings.Sc=17;          % Schmidt number
settings.Gamma=350;      % Dimensionless gravitational force (Richardson number)
settings.Dr=2.12;        % Dimensionless Rotational Coefficient
settings.Q=0.5;          % Flow rate setting

%% ------ Mesh Initialization ---------------------------------------------
N_mesh=175; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,4,bc_type.none,tran_type.none);
% mesh=cfd6(N_mesh,4,bc_type.none);

%% ------ Solve for basic state--------------------------------------------
% Solving for base flow
settings.disp_iter=true;

% prof_obj_G0=osp_Hbase_const_G0(cheb,settings);
% prof_obj_G0=prof_obj_G0.solve();
% 
% U0_G0=prof_obj_G0.xout(2:N_mesh+1);
% H0_G0=prof_obj_G0.xout(N_mesh+2:2*N_mesh+1);
% 
% settings.Gamma=prof_obj_G0.xout(1);

prof_obj=osp_Hbase_const(mesh,settings);
prof_obj=prof_obj.solve();

U0=prof_obj.xout(2:N_mesh+1);
H0=prof_obj.xout(N_mesh+2:2*N_mesh+1);

%% Plotting
figure(3);
subplot(2,1,1);
% hold on;
plot(mesh.col_pt,U0);
% plot(cheb.col_pt,U0_G0);
% hold off;

subplot(2,1,2);
% hold on;
plot(mesh.col_pt,H0);
% plot(cheb.col_pt,H0_G0);
% hold off;


