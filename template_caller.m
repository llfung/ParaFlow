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
settings.Re=0.126;       % Reynolds number based on half width and swimming speed
settings.Gamma=180;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=0.5;
settings.Q=-4;          % Flow rate setting

%% ------ Mesh Initialization ---------------------------------------------
N_mesh=175; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,4,bc_type.none,tran_type.none);
% mesh=cfd6(N_mesh,4,bc_type.none);

%% ------ Solve for basic state--------------------------------------------
% Solving for base flow
settings.disp_iter=true;

prof_obj=NSTrans_VDlib(mesh,settings);
prof_obj=prof_obj.solve();

U0=prof_obj.prof(2:N_mesh+1);
N0=prof_obj.prof(N_mesh+2:2*N_mesh+1);

%% Plotting
figure(3);
subplot(2,1,1);
% hold on;
plot(mesh.col_pt,U0);
% plot(cheb.col_pt,U0_G0);
% hold off;

subplot(2,1,2);
% hold on;
plot(mesh.col_pt,N0);
% plot(cheb.col_pt,H0_G0);
% hold off;


