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
settings.Gamma=0.1;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=1;
settings.Q=-0.5;          % Flow rate setting
%% ------ Mesh Initialization ---------------------------------------------
N_mesh=101; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,2,bc_type.none,tran_type.lin);
% mesh=cfd6(N_mesh,4,bc_type.none);

%% ------ Profile Initialization ------------------------------------------
G0=0;
U0=(1-mesh.pts.^2)*3/4*settings.Q;
%U0=cos(3*pi*mesh.pts/2)*3*pi/2*settings.Q;
N0=ones(N_mesh,1);

x0=[G0;U0;N0];
%% ------ Solve for basic state--------------------------------------------
% Solving for base flow
settings.disp_iter=true;

prof_obj=NSTransH_VDlib_sym(mesh,settings,x0);
prof_obj=prof_obj.solve();

U0=prof_obj.prof(2:N_mesh+1);
N0=prof_obj.prof(N_mesh+2:2*N_mesh+1);

%% Plotting
figure(1);
subplot(2,1,1);
% hold on;
plot(mesh.col_pt,U0);
ylabel('U(x)');
% plot(cheb.col_pt,U0_G0);
% hold off;

subplot(2,1,2);
% hold on;
plot(mesh.col_pt,N0);
xlabel('x');
ylabel('N(x)');
% plot(cheb.col_pt,H0_G0);
% hold off;


