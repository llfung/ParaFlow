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
settings.Gamma=4.5;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=0.5;      % Swimming Peclet Number
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
settings.limiter_value=1000; % Stop iteration if N(0)>1000
settings.limiter_index=2*N_mesh+1; % Index for N(0)


cont_obj_pre=NSTransH_VDconst_contRi_G0_sym(mesh,settings,0.2,2000,1);
cont_obj_pre=cont_obj_pre.run_cont_forward();
settings.Gamma=cont_obj_pre.para_out(1,end);
settings.G0=0;
cont_obj=NSTransH_VDconst_contRi_G0_sym(mesh,settings,0.2,3500,1,cont_obj_pre.x_out(:,end));
cont_obj=cont_obj.run_cont_backward();

%% Some Post-Processing
x=mesh.col_pt;
Ux=cont_obj.x_out(2:(N_mesh+1),:);
Q=mesh.wint*Ux;
Ux=Ux-Q*(mesh.wint*ones(N_mesh,1));
Nx=exp(cont_obj.x_out(N_mesh+2:(2*N_mesh+1),:))/2;

%% U(r)/ N(r) plot
U0_select=-1:-1:-5;
col_order=[0:1/(length(U0_select)-1):1]'*[1 0 0]...
    +[1:-1/(length(U0_select)-1):0]'*[0 0 1];

Uxf=figure('Position',[100,320,450,275]);hold on;
colororder(Uxf,col_order);Uxa=gca;
set(Uxa,'FontSize',14);
Nxf=figure('Position',[400,320,450,275]);hold on;
colororder(Nxf,col_order);Nxa=gca;
set(Nxa,'FontSize',14);

leg=cell(1,length(U0_select));
for i=1:length(U0_select)
    ind=find(Ux(N_mesh,:)<U0_select(i),1,'first');
        plot(Uxa,[-x;x(end-1:-1:1)],[Ux(:,ind);Ux(end-1:-1:1,ind)],'-');
        plot(Nxa,[-x;x(end-1:-1:1)],[Nx(:,ind);Nx(end-1:-1:1,ind)],'-');
    leg{i}=['Ri=' num2str(cont_obj.para_out(ind),4)];
end
figure(Uxf);hold off;
xlabel('x');ylabel('u(x)');
% legend(leg,'location','south');
figure(Nxf);hold off;
% legend(leg,'location','south');
xlabel('x');ylabel('n(x)');
% set(Nxa,'YScale','log');

% saveas(Uxf,'Ux.fig','fig');
% saveas(Uxf,'Ux.eps','epsc');
% saveas(Nxf,'Nx.fig','fig');
% saveas(Nxf,'Nx.eps','epsc');

%% U(0)-Ri Plot
Bi_f=figure('Position',[700,320,350,320]);
hold on;Bifur_a=gca;
set(Bifur_a,'FontSize',14);
% plot(cont_obj.para_out,cont_obj.x_out(1+floor((N_mesh+1)/2),:));
plot(Bifur_a,cont_obj.para_out,Ux(N_mesh,:),'b-');
plot(Bifur_a,[0 pi^2/(1.1*settings.Re/settings.Pe_s)],[0 0],'b-');
plot(Bifur_a,[pi^2/(1.1*settings.Re/settings.Pe_s) 200],[0 0],'b--');

for i=1:length(U0_select)
    ind=find(Ux(N_mesh,:)<U0_select(i),1,'first');
    plot(cont_obj.para_out(ind),Ux(N_mesh,ind),'o','Color',col_order(i,:));
end
hold off;

% legend('Amplitude u(0)','Ri Re \xi n_0/Pe_s= \pi^2 (1+(u(0) \xi / Pe_s)^2/12)','location','northwest','FontSize',11);
xlabel('$$Ri$$','Interpreter','latex','FontSize',14);
ylabel('$$u(0)$$','Interpreter','latex','FontSize',14)
axis([0 150 -6 5]);

% saveas(Bifur_f,'sedi_bi.fig','fig');
% saveas(Bifur_f,'sedi_bi.eps','epsc');
%% U(0)-Ri Plot
Bifur_f=figure('Position',[700,320,350,320]);
hold on;Bifur_a=gca;
set(Bifur_a,'FontSize',14);
% plot(cont_obj.para_out,cont_obj.x_out(1+floor((N_mesh+1)/2),:));
plot(Bifur_a,cont_obj.para_out *0.93* settings.Re/settings.Pe_s,Ux(N_mesh,:) /settings.Pe_s,'b-');
plot(Bifur_a,pi^2+[-10:0.01:10].^2*pi^2/12,[-10:0.01:10],'r:');
plot(Bifur_a,[0 pi^2],[0 0],'b-');
plot(Bifur_a,[pi^2 100],[0 0],'b--');

for i=1:length(U0_select)
    ind=find(Ux(N_mesh,:)<U0_select(i),1,'first');
    plot(cont_obj.para_out(ind)*0.93*settings.Re/settings.Pe_s,Ux(N_mesh,ind)/settings.Pe_s,'o','Color',col_order(i,:));
end
hold off;

legend('Amplitude u(0)','Ri Re \xi n_0/Pe_s= \pi^2 (1+(u(0) \xi / Pe_s)^2/12)','location','northwest','FontSize',11);
xlabel('$$Ri Re \xi n_0/Pe_s$$','Interpreter','latex','FontSize',14);
ylabel('$$u(0) \xi / Pe_s$$','Interpreter','latex','FontSize',14)
axis([0 20 -4 5]);

% saveas(Bifur_f,'sedi_bi.fig','fig');
% saveas(Bifur_f,'sedi_bi.eps','epsc');