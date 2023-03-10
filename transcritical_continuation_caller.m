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
settings.Gamma=4;      % Dimensionless gravitational force (Richardson number)
settings.Pe_s=1;      % Swimming Peclet Number
%settings.Q=0.01;          % Flow rate setting
settings.G0=-0.01;

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


cont_obj_pre=NSTransH_VDconst_contRi_G0_radial(mesh,settings,0.1,1000,1);
cont_obj_pre=cont_obj_pre.run_cont_forward();
settings.Gamma=cont_obj_pre.para_out(1,end);
settings.G0=0;
cont_obj=NSTransH_VDconst_contRi_G0_radial(mesh,settings,0.2,1000,1,cont_obj_pre.x_out(:,end));
cont_obj=cont_obj.run_cont_backward();

%% Some Post-Processing
wint=mesh.wint.*(mesh.col_pt');
r=mesh.col_pt;
Ur=cont_obj.x_out(2:(N_mesh+1),:);
Q=wint*Ur;
Ur=Ur-Q*(wint*ones(N_mesh,1));
Nr=exp(cont_obj.x_out(N_mesh+2:(2*N_mesh+1),:));
ind0=find(Ur(N_mesh,:)<0,1,'first');

%% U(r)/ N(r) plot
U0_select=0.5:-0.5:-3.5;
col_order=[0:1/(length(U0_select)-1):1]'*[1 0 0]...
    +[1:-1/(length(U0_select)-1):0]'*[0 0 1];

Urf=figure('Position',[100,320,450,275]);hold on;
colororder(Urf,col_order);Ura=gca;
set(Ura,'FontSize',14);
Nrf=figure('Position',[400,320,450,275]);hold on;
colororder(Nrf,col_order);Nra=gca;
set(Nra,'FontSize',14);

leg=cell(1,length(U0_select));
for i=1:length(U0_select)
    ind=find(Ur(N_mesh,:)<U0_select(i),1,'first');
    if Ur(N_mesh,ind-1)<0
        plot(Ura,r,Ur(:,ind),'--');
        plot(Nra,r,Nr(:,ind),'--');
    else
        plot(Ura,r,Ur(:,ind),'-');
        plot(Nra,r,Nr(:,ind),'-');
    end
    leg{i}=['Ri=' num2str(cont_obj.para_out(ind),4)];
end
figure(Urf);hold off;
xlabel('r');ylabel('u(r)');
%legend(leg,'location','southeast');
figure(Nrf);hold off;
%legend(leg,'location','northeast');
xlabel('r');ylabel('n(r)');
set(Nra,'YScale','log');

% saveas(Urf,'Ur.fig','fig');
% saveas(Urf,'Ur.eps','epsc');
% saveas(Nrf,'Nr.fig','fig');
% saveas(Nrf,'Nr.eps','epsc');

%% U(0)-Ri Plot
Bifur_f=figure('Position',[700,320,819,518]);
hold on;Bifur_a=gca;
set(Bifur_a,'FontSize',24);
% plot(Bifur_a,cont_obj.para_out,cont_obj.x_out(1+floor((N_mesh+1)/2),:));
plot(Bifur_a,cont_obj.para_out(1:ind0-1) *0.93* settings.Re/settings.Pe_s,Ur(N_mesh,1:ind0-1)*0.93 /settings.Pe_s,'b-','LineWidth',2);
plot(Bifur_a,[8 8],[-20 20],'r:','LineWidth',2);
plot(Bifur_a,cont_obj.para_out(ind0:end) *0.93* settings.Re/settings.Pe_s,Ur(N_mesh,ind0:end)*0.93 /settings.Pe_s,'b--','LineWidth',2);
%plot(Bifur_a,pi^2+[-10:0.01:10].^2*pi^2/12,[-10:0.01:10],'r:');
plot(Bifur_a,[0 14.68197064212389],[0 0],'b-','LineWidth',2);
plot(Bifur_a,[14.68197064212389 100],[0 0],'b--','LineWidth',2);

for i=1:length(U0_select)
    ind=find(Ur(N_mesh,:)<U0_select(i),1,'first');
    plot(cont_obj.para_out(ind)*0.93* settings.Re/settings.Pe_s,Ur(N_mesh,ind)*0.93/settings.Pe_s,'o','Color',col_order(i,:));
end
hold off;

legend('$$u(0)$$','$$3 Ri Re \alpha_0 n_0/ V_s =8$$','location','northwest','interpreter','latex','FontSize',18);
xlabel('$$3 Ri Re \alpha_0 n_0/ V_s$$','Interpreter','latex','FontSize',24);
ylabel('$$3 u(0) \alpha_0/  V_s$$','Interpreter','latex','FontSize',24)
axis([0 20 -8 2]);

% saveas(Bifur_f,'sedi_trans.fig','fig');
% saveas(Bifur_f,'sedi_trans.eps','epsc');
