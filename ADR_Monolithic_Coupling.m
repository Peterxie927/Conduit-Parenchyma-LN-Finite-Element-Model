clear all
close all
clc

%% Load Conduit Network into Workspace
tic
K_dimensionless =  [2.5E-08 2.5E-08];
D = 1E-10;
load("Conduit_microscale_highres_el_0_35.mat")
% load("Conduit_model_network_v6_Frontal_Delaunay.mat")
% load("Conduit_model_network_v5_Delaunay.mat");
% load("Conduit_model_network_v6_Frontal_Delaunay.mat")
width_x = 1E-04; % Width of domain (in x-direction) [m]
width_y = 1E-04; % Width of domain (in y-direction) [m]
NL = NL(:,1:2)*(1E-04);
EL = EL(:,1:3);
NoN = size(NL,1);
Ly = width_y/(round(sqrt(NoN))); p27 = 1; p10 = 1; p11 = 1;

for p9 = 1:size(Line_list,1)
    if ((NL(Line_list(p9,1),1)==0) && (NL(Line_list(p9,2),1) ==0)) || ((NL(Line_list(p9,1),1)==width_x) && (NL(Line_list(p9,2),1) ==width_x)) || ((NL(Line_list(p9,1),2)==0) && (NL(Line_list(p9,2),2) ==0)) || ((NL(Line_list(p9,1),2)==width_y) && (NL(Line_list(p9,2),2) == width_y))
        % ((ismember(Line_list(1),Left_indices)&&ismember(Line_list(2),Left_indices))||(ismember(Line_list(1),Right_indices)&&ismember(Line_list(2),Right_indices))||(ismember(Line_list(1),Top_indices)&&ismember(Line_list(2),Top_indices))||(ismember(Line_list(1),Bottom_indices)&&ismember(Line_list(2),Bottom_indices)))
    else
        Conduit_Surf(p27,:) = [Line_list(p9,1) Line_list(p9,2)];
        conduit_coordinates1(p27,:) = [NL(Line_list(p9,1),1)];
        conduit_coordinates2(p27,:) = [NL(Line_list(p9,1),2)];
        conduit_coordinates3(p27,:) = [NL(Line_list(p9,2),1)];
        conduit_coordinates4(p27,:) = [NL(Line_list(p9,2),2)];
        conduit_indices(p27) = p9;
        Line_list2(p27,:) = [Line_list(p9,1) Line_list(p9,2)];
        % bound_node1 are given source BCs, bound_node2 are given sink BCs.
        if (((NL(Line_list(p9,1),1)==width_x) && ((NL(Line_list(p9,1),2)==width_y))))
            bound_node1(p10) = [Line_list(p9,1)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==(10^(-4)*0.3))))
            bound_node2(p11) = [Line_list(p9,2)];
            p11 = p11+1;
        elseif (((NL(Line_list(p9,2),1)==0.3*10^(-4))) && ((NL(Line_list(p9,2),2)==0)))
            bound_node2(p11) = [Line_list(p9,2)];
            p11 = p11+1;
        end
        p27 = p27+1;
    end
end

[Pressure_sol] = Fluid_transport_solver(NL,EL,NoN,Ly,Line_list,width_x,width_y,Conduit_Surf, bound_node1, bound_node2,conduit_coordinates1,conduit_coordinates2,conduit_coordinates3,conduit_coordinates4,Line_list2);
[Ux Uy] = gradient_solver(Pressure_sol,NL,EL,NoN,Ly,width_x,width_y);

Pe = [(Ux.*width_x)./D (Uy.*width_y)./D];
figure (1)
subplot(2,2,1)
patch('Faces',EL(:,1:3),'Vertices',NL,'CData',Pressure_sol,'FaceColor','interp');
h=colorbar;
axis([0 width_x 0 width_y]);
title(h,'[ ]');
xlabel('x [m]','FontSize',16,'interpreter','Latex');
ylabel('y [m]','FontSize',16,'interpreter','Latex');
set(gcf,'color','w'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
% title('Pressure Plot from Seepage flow (S/S)')
legend('Pressure Field')
hold on
for i = 1:size(conduit_coordinates1,1)
    plot([conduit_coordinates1(i) conduit_coordinates3(i)],[conduit_coordinates2(i) conduit_coordinates4(i)],'LineWidth',1.8,'color','k')
    hold on
end
set(gca,'XTick',[]); %which will get rid of all the markings for the y axis
set(gca,'YTick',[]); %which will get rid of all the markings for the y axis

subplot(2,2,3)
patch('Faces',EL(:,1:3),'Vertices',NL,'CData',Ux,'FaceColor','flat');
hold on
h=colorbar;
title(h,'[ ]');
% xlabel('x [m]','FontSize',16,'interpreter','Latex');
% ylabel('y [m]','FontSize',16,'interpreter','Latex');
axis([0 width_x 0 width_y]); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gcf,'color','w')
title('UX (S/S)')
for i = 1:size(conduit_coordinates1,1)
    plot([conduit_coordinates1(i) conduit_coordinates3(i)],[conduit_coordinates2(i) conduit_coordinates4(i)],'LineWidth',1.8,'color','k')
    hold on
end
subplot(2,2,4)
patch('Faces',EL(:,1:3),'Vertices',NL,'CData',Uy,'FaceColor','flat');
h=colorbar;
title(h,'[ ]');
% xlabel('x [m]','FontSize',16,'interpreter','Latex');
% ylabel('y [m]','FontSize',16,'interpreter','Latex');
set(gcf,'color','w'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
title('UY (S/S)')
hold on
axis([0 width_x 0 width_y]);
for i = 1:size(conduit_coordinates1,1)
    plot([conduit_coordinates1(i) conduit_coordinates3(i)],[conduit_coordinates2(i) conduit_coordinates4(i)],'LineWidth',1.8,'color','k')
    hold on
end

[C4 C_ot dt total_time] = ADR_Transport_Solver(NL,EL,NoN,Ly,Pe,Line_list,width_x,width_y,Conduit_Surf, bound_node1, bound_node2,conduit_coordinates1,conduit_coordinates2,conduit_coordinates3,conduit_coordinates4,Line_list2);


%% Mesh Interpolation

[A,a,b,c] = ShapeFunctionCoeff(NL,EL);
i2 = 0;
total_time = 1E-02;
dt = 1E-05;
no_t_steps = round(total_time/dt); % Number of time-steps [ ]

% [X,Y] = meshgrid(linspace(0,width_x,101),linspace(0,width_y,101));
% x = linspace(0.5,width_x-0.05,5); y = linspace(0,width_y,101);
x = [0.5 0.7 0.9]*(1E-04); y = [0.4 0.6 0.8]*(1E-04);
P_array = zeros(no_t_steps/10,length(x),length(y));
tir = 0;

for ti = 1:10:no_t_steps
    tir = tir+1;
    for xi = 1:length(x)
        for yi = 1:length(y)
            P = [x(xi) y(yi)];
            for i=1:size(EL,1)
                P1 = NL(EL(i,1),:); P2 = NL(EL(i,2),:); P3 = NL(EL(i,3),:);
                % P12 = P1-P2; P23 = P2-P3; P31 = P3-P1;
                % t = (sign(det([P31;P23]))*sign(det([P3-P;P23])) >= 0) & (sign(det([P12;P31]))*sign(det([P1-P;P31])) >= 0) & (sign(det([P23;P12]))*sign(det([P2-P;P12])) >= 0);
                Ptri=[P1;P2;P3];
                if (P(1)<=max(Ptri(:,1)))&&(P(1)>=min(Ptri(:,1)))&&(P(2)<=max(Ptri(:,2)))&&(P(2)>=min(Ptri(:,2)))
                    i2 = i2+1;
                    EL_s(i2) = i;
                end
            end
            El_s = EL_s(1);
            % add the linear combination of nodal values weighted by shape functions
            % w_SF1 = (1/(2*A(El_s)))*(a(El_s,1)+b(El_s,1)*(P(1)-NL(EL(El_s,1),1))+c(El_s,1)*(P(2)-NL(EL(El_s,1),2)))*Pressure_sol(EL(El_s,1));
            % w_SF2 = (1/(2*A(El_s)))*(a(El_s,2)+b(El_s,2)*(P(1)-NL(EL(El_s,2),1))+c(El_s,2)*(P(2)-NL(EL(El_s,2),2)))*Pressure_sol(EL(El_s,2));
            % w_SF3 = (1/(2*A(El_s)))*(a(El_s,3)+b(El_s,3)*(P(1)-NL(EL(El_s,3),1))+c(El_s,3)*(P(2)-NL(EL(El_s,3),2)))*Pressure_sol(EL(El_s,3));
            w_SF1 = (1/(2*A(El_s)))*(a(El_s,1)+b(El_s,1)*(P(1))+c(El_s,1)*(P(2)))*C_ot(EL(El_s,1),ti);
            w_SF2 = (1/(2*A(El_s)))*(a(El_s,2)+b(El_s,2)*(P(1))+c(El_s,2)*(P(2)))*C_ot(EL(El_s,2),ti);
            w_SF3 = (1/(2*A(El_s)))*(a(El_s,3)+b(El_s,3)*(P(1))+c(El_s,3)*(P(2)))*C_ot(EL(El_s,3),ti);
    
            P_array(tir,xi,yi) = (w_SF1+w_SF2+w_SF3);
            i2 = 0; EL_s = 0;
        end
    end
end

save('Conduit_profiles_0_35.mat')