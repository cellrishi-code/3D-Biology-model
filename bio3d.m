%% ============================================================
%  BIOLOGICAL CELL SIGNALING — 2D + 3D SIMULATION
%  Figure 1 (2D): HH Action Potential, Gates, Conductances,
%                 Ca2+ Oscillations, Receptor Binding, Phase Portrait
%  Figure 2 (3D): Phase Space Trajectory, Ca2+ Wave Grid,
%                 Axon Propagation Surface, Receptor Landscape,
%                 Conductance Torus, Full Cell State Cloud
%% ============================================================
clear; clc; close all;

%% ---- Settings ----------------------------------------------
dt=0.01; tEnd=100; t=0:dt:tEnd; N=length(t);

%% ============================================================
%  MODULE 1 — HODGKIN-HUXLEY
%% ============================================================
Cm=1; gNa=120; gK=36; gL=0.3;
ENa=50; EK=-77; EL=-54.4; Vrest=-65;

I_ext=zeros(1,N);
I_ext(t>=10 & t<=11)=10;
I_ext(t>=40 & t<=41)=10;
I_ext(t>=60 & t<=60.5)=14;
I_ext(t>=75 & t<=76)=12;

alpha_m=@(V) 0.1*(V+40)./(1-exp(-(V+40)/10));
beta_m =@(V) 4*exp(-(V+65)/18);
alpha_h=@(V) 0.07*exp(-(V+65)/20);
beta_h =@(V) 1./(1+exp(-(V+35)/10));
alpha_n=@(V) 0.01*(V+55)./(1-exp(-(V+55)/10));
beta_n =@(V) 0.125*exp(-(V+65)/80);

V_hh=zeros(1,N); V_hh(1)=Vrest;
m=zeros(1,N); m(1)=alpha_m(Vrest)/(alpha_m(Vrest)+beta_m(Vrest));
h=zeros(1,N); h(1)=alpha_h(Vrest)/(alpha_h(Vrest)+beta_h(Vrest));
n=zeros(1,N); n(1)=alpha_n(Vrest)/(alpha_n(Vrest)+beta_n(Vrest));
gNa_t=zeros(1,N); gK_t=zeros(1,N);

for i=1:N-1
    Vm=V_hh(i);
    INa=gNa*m(i)^3*h(i)*(Vm-ENa);
    IK =gK*n(i)^4*(Vm-EK);
    IL =gL*(Vm-EL);
    gNa_t(i)=gNa*m(i)^3*h(i);
    gK_t(i) =gK*n(i)^4;
    V_hh(i+1)=Vm+dt/Cm*(I_ext(i)-INa-IK-IL);
    m(i+1)=m(i)+dt*(alpha_m(Vm)*(1-m(i))-beta_m(Vm)*m(i));
    h(i+1)=h(i)+dt*(alpha_h(Vm)*(1-h(i))-beta_h(Vm)*h(i));
    n(i+1)=n(i)+dt*(alpha_n(Vm)*(1-n(i))-beta_n(Vm)*n(i));
end
gNa_t(end)=gNa_t(end-1); gK_t(end)=gK_t(end-1);

%% ============================================================
%  MODULE 2 — CALCIUM OSCILLATIONS (Li-Rinzel)
%% ============================================================
c0=2; c1=0.185; v1=6; v2=0.11; v3=0.9;
d1=0.13; d2=1.049; d3=0.9434; d5=0.08234; a2=0.2; IP3=0.5;

Ca=zeros(1,N); Ca(1)=0.1;
Ce=zeros(1,N); Ce(1)=(c0-Ca(1))/c1;
h_ca=zeros(1,N); h_ca(1)=0.85;

for i=1:N-1
    Cai=Ca(i); Cei=Ce(i); hi=h_ca(i);
    m_inf=IP3/(IP3+d1); n_inf=Cai/(Cai+d5);
    J_chan=v1*(m_inf*n_inf*hi)^3*(Cei-Cai);
    J_leak=v2*(Cei-Cai);
    J_pump=v3*Cai^2/(Cai^2+0.3^2);
    h_inf=d2*(IP3+d1)/((d2*(IP3+d1))+a2*(IP3+d3));
    tau_h=1/(a2*(IP3+d3)+a2*Cai);
    Ca(i+1)  =max(Cai+dt*(J_chan+J_leak-J_pump),0.001);
    Ce(i+1)  =max(Cei+dt*(-J_chan-J_leak+J_pump)/c1,0.001);
    h_ca(i+1)=hi+dt*(h_inf-hi)/tau_h;
end

%% ============================================================
%  MODULE 3 — RECEPTOR-LIGAND BINDING
%% ============================================================
kon=0.5; koff=0.05; R_total=0.5; L_total=1.0;
L_conc=zeros(1,N); L_conc(t>=20)=L_total; L_conc(t>=70)=0;
RL=zeros(1,N); R=zeros(1,N); R(1)=R_total;
for i=1:N-1
    dRL=kon*R(i)*L_conc(i)-koff*RL(i);
    RL(i+1)=max(0,min(RL(i)+dt*dRL,R_total));
    R(i+1)=R_total-RL(i+1);
end
occ=RL/R_total;

%% ============================================================
%  MODULE 4 — NOISE
%% ============================================================
rng(42); V_noise=zeros(1,N); sigma=1.5; tau_n=2.0;
for i=1:N-1
    V_noise(i+1)=V_noise(i)*exp(-dt/tau_n)+sigma*sqrt(1-exp(-2*dt/tau_n))*randn();
end
V_noisy=V_hh+V_noise;

%% ============================================================
%  MODULE 5 — SPATIAL CALCIUM WAVE (40x40 cell grid)
%% ============================================================
GX=40; GY=40; D_ca=0.05;
tskip=50; t_spa=t(1:tskip:end); Nt_spa=length(t_spa);
Ca_grid=zeros(GX,GY,Nt_spa);
Ca2d=0.1*ones(GX,GY);
cx=GX/2; cy=GY/2;
[gxi,gyi]=ndgrid(1:GX,1:GY);
dist2=((gxi-cx).^2+(gyi-cy).^2);

for k=1:Nt_spa
    Ca_local=Ca(min(round(t_spa(k)/dt)+1,N));
    Ca2d=Ca_local*(0.6+0.4*exp(-dist2/100).*(1+0.3*sin(gxi/3)+0.3*cos(gyi/3)));
    Ca2d=max(Ca2d,0.001);
    lap=[Ca2d(2:end,:);Ca2d(end,:)]+[Ca2d(1,:);Ca2d(1:end-1,:)]+...
        [Ca2d(:,2:end),Ca2d(:,end)]+[Ca2d(:,1),Ca2d(:,1:end-1)]-4*Ca2d;
    Ca2d=Ca2d+D_ca*lap;
    Ca_grid(:,:,k)=max(Ca2d,0.001);
end

%% ============================================================
%  MODULE 6 — AXON PROPAGATION (1D cable, HH)
%% ============================================================
Nx=60; dx=10; r_a=0.1;
Nx_t=500; dt_ax=0.05; t_ax=0:dt_ax:(Nx_t-1)*dt_ax;
V_axon=Vrest*ones(Nx,Nx_t);
m_ax=m(1)*ones(Nx,Nx_t); h_ax=h(1)*ones(Nx,Nx_t); n_ax=n(1)*ones(Nx,Nx_t);

for i=1:Nx_t-1
    V=V_axon(:,i); ma=m_ax(:,i); ha=h_ax(:,i); na=n_ax(:,i);
    Iax=zeros(Nx,1);
    Iax(2:end-1)=(V(1:end-2)-2*V(2:end-1)+V(3:end))/(r_a*dx^2);
    Istim=zeros(Nx,1);
    if t_ax(i)>=2  && t_ax(i)<=2.5;  Istim(5)=12; end
    if t_ax(i)>=15 && t_ax(i)<=15.5; Istim(5)=12; end
    if t_ax(i)>=28 && t_ax(i)<=28.5; Istim(5)=12; end
    INa_ax=gNa*ma.^3.*ha.*(V-ENa);
    IK_ax =gK*na.^4.*(V-EK);
    IL_ax =gL*(V-EL);
    V_axon(:,i+1)=V+dt_ax/Cm*(Iax+Istim-INa_ax-IK_ax-IL_ax);
    m_ax(:,i+1)=ma+dt_ax*(alpha_m(V).*(1-ma)-beta_m(V).*ma);
    h_ax(:,i+1)=ha+dt_ax*(alpha_h(V).*(1-ha)-beta_h(V).*ha);
    n_ax(:,i+1)=na+dt_ax*(alpha_n(V).*(1-na)-beta_n(V).*na);
end

%% ============================================================
%  STYLING HELPERS
%% ============================================================
cbg=[0.07 0.07 0.12]; cax=[0.15 0.15 0.22]; cg=[0.25 0.25 0.35];
C={[0.2 0.8 1],[1 0.4 0.4],[0.4 1 0.4],[1 0.85 0.2],[0.8 0.4 1],[1 0.6 0.2]};

function style2d(ax,cax,cg)
    set(ax,'Color',cax,'XColor',[0.8 0.8 0.9],'YColor',[0.8 0.8 0.9],...
        'GridColor',cg,'GridAlpha',0.4,'Box','off','FontSize',9);
    grid(ax,'on');
end

function style3d(ax,cax,cg)
    set(ax,'Color',cax,'XColor',[0.8 0.8 0.9],'YColor',[0.8 0.8 0.9],...
        'ZColor',[0.8 0.8 0.9],'GridColor',cg,'GridAlpha',0.4,'Box','off','FontSize',9);
    grid(ax,'on');
end

%% ============================================================
%  FIGURE 1 — 2D PANELS
%% ============================================================
fig1=figure('Name','2D Cell Signals','Color',cbg,'Position',[20 20 1400 900]);

% ① Action Potential
ax1=subplot(4,2,[1 2]); style2d(ax1,cax,cg); hold on;
plot(t,V_noisy,'Color',[0.5 0.5 0.6 0.5],'LineWidth',0.5);
plot(t,V_hh,'Color',C{1},'LineWidth',1.8);
fill([t fliplr(t)],[I_ext/max(I_ext)*15-84, -84*ones(1,N)],C{4},'FaceAlpha',0.3,'EdgeColor','none');
ylabel('V_m (mV)','Color',[0.9 0.9 1]); ylim([-90 60]); xlim([0 tEnd]);
title('① Action Potential (Hodgkin-Huxley)','Color','w','FontSize',12,'FontWeight','bold');
legend({'Noisy','V_m','Stimulus'},'TextColor',[0.9 0.9 1],'Color',cax,'EdgeColor',cg);

% ② Gate variables
ax2=subplot(4,2,3); style2d(ax2,cax,cg); hold on;
plot(t,m,'Color',C{1},'LineWidth',1.5);
plot(t,h,'Color',C{2},'LineWidth',1.5);
plot(t,n,'Color',C{3},'LineWidth',1.5);
ylabel('Probability','Color',[0.9 0.9 1]); ylim([0 1]); xlim([0 tEnd]);
title('② Ion Channel Gate Variables','Color','w','FontSize',11,'FontWeight','bold');
legend({'m (Na act.)','h (Na inact.)','n (K act.)'},'TextColor',[0.9 0.9 1],'Color',cax,'EdgeColor',cg,'FontSize',8);

% ③ Conductances
ax3=subplot(4,2,4); hold on; grid on;
set(ax3,'Color',cax,'GridColor',cg,'GridAlpha',0.4,'Box','off','FontSize',9);
yyaxis left;  plot(t,gNa_t,'Color',C{1},'LineWidth',1.5); ylabel('g_{Na}','Color',C{1}); ax3.YAxis(1).Color=C{1};
yyaxis right; plot(t,gK_t,'Color',C{3},'LineWidth',1.5);  ylabel('g_K','Color',C{3});  ax3.YAxis(2).Color=C{3};
ax3.XColor=[0.8 0.8 0.9]; xlim([0 tEnd]);
title('③ Ionic Conductances','Color','w','FontSize',11,'FontWeight','bold');

% ④ Calcium
ax4=subplot(4,2,[5 6]); style2d(ax4,cax,cg); hold on;
area(t,Ca,'FaceColor',C{2},'FaceAlpha',0.2,'EdgeColor','none');
plot(t,Ca,'Color',C{2},'LineWidth',2);
plot(t,Ce*0.1,'--','Color',[1 0.7 0.5],'LineWidth',1.2);
plot(t,h_ca,'Color',C{5},'LineWidth',1.2);
ylabel('[Ca²⁺] µM','Color',[0.9 0.9 1]); xlim([0 tEnd]);
title('④ Ca²⁺ Oscillations (Li-Rinzel)','Color','w','FontSize',11,'FontWeight','bold');
legend({'area','[Ca]_{cyt}','[Ca]_{ER}×0.1','h_ca'},'TextColor',[0.9 0.9 1],'Color',cax,'EdgeColor',cg,'FontSize',8);

% ⑤ Receptor binding
ax5=subplot(4,2,7); style2d(ax5,cax,cg); hold on;
fill([t fliplr(t)],[occ zeros(1,N)],C{4},'FaceAlpha',0.35,'EdgeColor','none');
plot(t,occ,'Color',C{4},'LineWidth',2);
plot(t,L_conc/L_total,'--','Color',[0.7 0.7 0.8],'LineWidth',1);
xline(20,'--','Color',C{2},'Alpha',0.8); text(21,0.9,'ON','Color',C{2},'FontSize',8);
xline(70,'--','Color',C{2},'Alpha',0.8); text(71,0.9,'OFF','Color',C{2},'FontSize',8);
ylabel('Fractional Occupancy','Color',[0.9 0.9 1]); xlabel('Time (ms)','Color',[0.9 0.9 1]);
title('⑤ Receptor-Ligand Binding','Color','w','FontSize',11,'FontWeight','bold');
ylim([0 1.1]); xlim([0 tEnd]);

% ⑥ Phase portrait
ax6=subplot(4,2,8); style2d(ax6,cax,cg); hold on;
cmap_pp=parula(N);
for ii=1:50:N-50; plot(V_hh(ii:ii+49),m(ii:ii+49),'Color',cmap_pp(ii,:),'LineWidth',1.2); end
scatter(V_hh(1),m(1),80,'g','filled','MarkerEdgeColor','w');
scatter(V_hh(end),m(end),80,'r','filled','MarkerEdgeColor','w');
xlabel('V_m (mV)','Color',[0.9 0.9 1]); ylabel('m gate','Color',[0.9 0.9 1]);
title('⑥ Phase Portrait (V_m vs m)','Color','w','FontSize',11,'FontWeight','bold');

sgtitle('BIOLOGICAL CELL SIGNAL SIMULATION — 2D','Color','w','FontSize',14,'FontWeight','bold');
set(fig1,'Color',cbg);

%% ============================================================
%  FIGURE 2 — 3D PANELS
%% ============================================================
fig2=figure('Name','3D Cell Signals','Color',cbg,'Position',[60 20 1440 950]);

%% ⑦ 3D Phase Space (Vm, m, n)
sp1=subplot(2,3,1); style3d(sp1,cax,cg); hold on;
skip3=5; idx3=1:skip3:N;
cmap3=cool(length(idx3));
for ii=1:length(idx3)-1
    i1=idx3(ii); i2=idx3(ii+1);
    plot3(V_hh(i1:i2),m(i1:i2),n(i1:i2),'Color',cmap3(ii,:),'LineWidth',1.5);
end
sp_i=find(V_hh>20);
if ~isempty(sp_i)
    scatter3(V_hh(sp_i),m(sp_i),n(sp_i),10,[1 0.5 0.2],'filled','MarkerFaceAlpha',0.6);
end
scatter3(V_hh(1),m(1),n(1),120,'g','filled','MarkerEdgeColor','w');
scatter3(V_hh(end),m(end),n(end),120,'r','filled','MarkerEdgeColor','w');
xlabel('V_m (mV)','Color',[0.9 0.9 1]); ylabel('m','Color',[0.9 0.9 1]); zlabel('n','Color',[0.9 0.9 1]);
title({'⑦ 3D Phase Space','Trajectory (V_m, m, n)'},'Color','w','FontSize',11,'FontWeight','bold');
colormap(sp1,cool); cb1=colorbar('Color',[0.8 0.8 0.9]);
ylabel(cb1,'Time \rightarrow','Color',[0.8 0.8 0.9],'FontSize',8);
view([-35 25]);

%% ⑧ 3D Calcium Wave Grid (stacked time frames)
sp2=subplot(2,3,2); style3d(sp2,cax,cg); hold on;
frames=[1, round(Nt_spa*0.25), round(Nt_spa*0.6), Nt_spa];
fcolors={[0.2 0.4 0.8],[0.2 0.7 0.9],[0.9 0.5 0.2],[1.0 0.2 0.2]};
falphas=[0.6 0.65 0.7 0.75];
[Xg,Yg]=meshgrid(1:GX,1:GY);
z_off=[0 0.6 1.2 1.8];
for fi=1:4
    Zg=Ca_grid(:,:,frames(fi))';
    surf(Xg,Yg,Zg+z_off(fi),'FaceAlpha',falphas(fi),'FaceColor',fcolors{fi},'EdgeAlpha',0.07);
    text(43,1,z_off(fi)+0.1,sprintf('t=%.0fms',t_spa(frames(fi))),'Color','w','FontSize',7,'FontWeight','bold');
end
xlabel('X (µm)','Color',[0.9 0.9 1]); ylabel('Y (µm)','Color',[0.9 0.9 1]); zlabel('[Ca²⁺]','Color',[0.9 0.9 1]);
title({'⑧ Ca²⁺ Wave','Propagation (Cell Grid)'},'Color','w','FontSize',11,'FontWeight','bold');
view([-40 30]);

%% ⑨ 3D Axon Propagation Surface
sp3=subplot(2,3,3); style3d(sp3,cax,cg); hold on;
t_ds=1:10:Nx_t; x_ds=1:2:Nx;
[Tx3,Xx3]=meshgrid(t_ax(t_ds),(x_ds-1)*dx);
surf(Tx3,Xx3,V_axon(x_ds,t_ds),'EdgeAlpha',0.04,'FaceAlpha',0.88);
xlabel('Time (ms)','Color',[0.9 0.9 1]); ylabel('Position (µm)','Color',[0.9 0.9 1]); zlabel('V_m (mV)','Color',[0.9 0.9 1]);
title({'⑨ 3D Axon Propagation','Spatiotemporal V_m'},'Color','w','FontSize',11,'FontWeight','bold');
colormap(sp3,jet); shading(sp3,'interp'); view([-50 30]);
cb3=colorbar('Color',[0.8 0.8 0.9]); ylabel(cb3,'V_m (mV)','Color',[0.8 0.8 0.9],'FontSize',8);

%% ⑩ 3D Receptor Occupancy Landscape
sp4=subplot(2,3,4); style3d(sp4,cax,cg); hold on;
kon_v=linspace(0.05,2,50); koff_v=linspace(0.01,0.5,50);
[KON,KOFF]=meshgrid(kon_v,koff_v);
occ_land=KON*L_total./(KON*L_total+KOFF);
surf(KON,KOFF,occ_land,'EdgeAlpha',0.08,'FaceAlpha',0.9);
our_occ=kon*L_total/(kon*L_total+koff);
scatter3(kon,koff,our_occ,160,'w','filled','MarkerEdgeColor',[1 0.85 0.2],'LineWidth',2);
text(kon+0.05,koff,our_occ+0.04,'Our cell','Color',[1 0.85 0.2],'FontSize',9,'FontWeight','bold');
xlabel('k_{on}','Color',[0.9 0.9 1]); ylabel('k_{off}','Color',[0.9 0.9 1]); zlabel('Occupancy','Color',[0.9 0.9 1]);
title({'⑩ Receptor Occupancy','Landscape (k_{on} vs k_{off})'},'Color','w','FontSize',11,'FontWeight','bold');
colormap(sp4,parula); shading(sp4,'interp'); view([-40 30]);

%% ⑪ 3D Conductance Torus
sp5=subplot(2,3,5); style3d(sp5,cax,cg); hold on;
theta=linspace(0,2*pi,80); phi=linspace(0,pi,80);
[TH,PH]=meshgrid(theta,phi);
gNa_norm=gNa_t/max(gNa_t);
t_idx=round(linspace(1,N,length(theta)));
gNa_s=gNa_norm(t_idx);
R_major=2; R_minor_s=0.3+0.7*gNa_s;
Xtor=(R_major + R_minor_s'.*cos(PH)).*cos(TH);
Ytor=(R_major + R_minor_s'.*cos(PH)).*sin(TH);
Ztor=R_minor_s'.*sin(PH);
Cdat=repmat(gNa_s',1,length(phi))';
surf(Xtor,Ytor,Ztor,Cdat,'EdgeAlpha',0.05,'FaceAlpha',0.88);
xlabel('X','Color',[0.9 0.9 1]); ylabel('Y','Color',[0.9 0.9 1]); zlabel('Z','Color',[0.9 0.9 1]);
title({'⑪ Conductance Torus','g_{Na} Modulates Shape'},'Color','w','FontSize',11,'FontWeight','bold');
colormap(sp5,turbo); shading(sp5,'interp'); axis(sp5,'equal'); view([30 25]);
cb5=colorbar('Color',[0.8 0.8 0.9]); ylabel(cb5,'g_{Na} (norm.)','Color',[0.8 0.8 0.9],'FontSize',8);

%% ⑫ 3D Cell State Cloud (Vm, Ca, Receptor Occ)
sp6=subplot(2,3,6); style3d(sp6,cax,cg); hold on;
Ca_full=interp1(t,Ca,t,'linear','extrap');
occ_full=interp1(t,occ,t,'linear','extrap');
skip_s=8; idx_s=1:skip_s:N;
Vs=V_hh(idx_s); Cs=Ca_full(idx_s); Os=occ_full(idx_s);
clr_t=linspace(0,1,length(idx_s));
scatter3(Vs,Cs,Os,18,clr_t,'filled','MarkerFaceAlpha',0.7);
plot3(Vs,Cs,Os,'Color',[1 1 1 0.12],'LineWidth',0.5);
rest_i=find(V_hh<-60 & Ca_full<0.15,5,'first');
if ~isempty(rest_i)
    scatter3(V_hh(rest_i(1)),Ca_full(rest_i(1)),occ_full(rest_i(1)),...
        140,'g','filled','MarkerEdgeColor','w','LineWidth',1.5);
    text(V_hh(rest_i(1))+1,Ca_full(rest_i(1)),occ_full(rest_i(1))+0.02,...
        'Rest','Color',[0.4 1 0.4],'FontSize',8,'FontWeight','bold');
end
pk_i=find(V_hh>30);
if ~isempty(pk_i)
    scatter3(V_hh(pk_i(1)),Ca_full(pk_i(1)),occ_full(pk_i(1)),...
        140,'r','filled','MarkerEdgeColor','w','LineWidth',1.5);
    text(V_hh(pk_i(1))+1,Ca_full(pk_i(1)),occ_full(pk_i(1))+0.02,...
        'Spike','Color',[1 0.4 0.4],'FontSize',8,'FontWeight','bold');
end
xlabel('V_m (mV)','Color',[0.9 0.9 1]); ylabel('[Ca²⁺] µM','Color',[0.9 0.9 1]); zlabel('Receptor Occ.','Color',[0.9 0.9 1]);
title({'⑫ Full Cell State Cloud','(V_m, [Ca²⁺], Occupancy)'},'Color','w','FontSize',11,'FontWeight','bold');
colormap(sp6,parula); view([-50 20]);
cb6=colorbar('Color',[0.8 0.8 0.9]); ylabel(cb6,'Time \rightarrow','Color',[0.8 0.8 0.9],'FontSize',8);

sgtitle('BIOLOGICAL CELL SIGNAL SIMULATION — 3D VIEWS','Color','w','FontSize',14,'FontWeight','bold');
set(fig2,'Color',cbg);

% Add lights to all 3D subplots
for ax_i=[sp1 sp2 sp3 sp4 sp5 sp6]
    axes(ax_i);
    light('Position',[1 1 2],'Style','infinite','Color',[0.8 0.8 1]);
    light('Position',[-1 -1 0],'Style','infinite','Color',[0.2 0.2 0.4]);
    lighting gouraud;
end

fprintf('\n=== Done! Fig 1 = 2D | Fig 2 = 3D ===\n');
fprintf('Peak Vm: %.1f mV | Spikes: %d\n', max(V_hh), sum(diff(V_hh>0)==1));
fprintf('Peak Ca: %.3f uM | Max Occ: %.1f%%\n', max(Ca), max(occ)*100);