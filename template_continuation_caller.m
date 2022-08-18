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
settings.Gamma=-10;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=1;      % Swimming Peclet Number
settings.Q=-.01;          % Flow rate setting

%% ------ Mesh Initialization ---------------------------------------------
N_mesh=175; %Number of mesh points (in r direction)

% Generating Chebyshev Collocation Points and Differential Matrices
mesh=chebyshev(N_mesh,4,bc_type.none,tran_type.lin);
% mesh=cfd6(N_mesh,2,bc_type.none);

%% ------ Continuation ----------------------------------------------------
settings.disp_iter=true;
settings.epsilon=1e-8;
% Limiter (limiting the value of N(0).
settings.limiter_value=1000; % Stop iteration if N(0)>1000
settings.limiter_index=2*N_mesh+1; % Index for N(0)

cont_obj=NSTransH_VDlib_contRi_sym(mesh,settings,.1,1000,1);
cont_obj=cont_obj.run_cont_forward();

%% Plotting
figure(4);hold on;
%plot(cont_obj.para_out,cont_obj.x_out(1+N_mesh+floor((N_mesh+1)/2),:));
plot(cont_obj.para_out,cont_obj.x_out(1+N_mesh,:));

figure(2);
subplot(2,1,1);
% hold on;
plot(mesh.col_pt,cont_obj.x_out(2:N_mesh+1,1:10:end));
% plot(cheb.col_pt,U0_G0);
% hold off;

subplot(2,1,2);
% hold on;
plot(mesh.col_pt,cont_obj.x_out(N_mesh+2:2*N_mesh+1,1:10:end));
% plot(cheb.col_pt,H0_G0);
% hold off;


