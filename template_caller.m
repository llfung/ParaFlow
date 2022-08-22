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
settings.Re=1;       % Reynolds number based on half width and swimming speed
settings.Gamma=1;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=1;
settings.Q=-0.01;          % Flow rate setting
%% ------ Mesh Initialization ---------------------------------------------
N_mesh=101; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,2,bc_type.none,tran_type.none);
% mesh=cfd6(N_mesh,4,bc_type.none);

%% ------ Profile Initialization ------------------------------------------
G0=0;
U0=(1-mesh.pts.^2)*3/2*settings.Q;
H0=ones(N_mesh,1);

x0=[G0;U0;H0];
%% ------ Solve for basic state--------------------------------------------
% Solving for base flow
settings.disp_iter=true;

prof_obj=NSTrans_VDlib(mesh,settings,x0);
prof_obj=prof_obj.solve();

U0=prof_obj.prof(2:N_mesh+1);
H0=prof_obj.prof(N_mesh+2:2*N_mesh+1);

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
plot(mesh.col_pt,exp(H0));
xlabel('x');
ylabel('N(x)');
% plot(cheb.col_pt,H0_G0);
% hold off;


