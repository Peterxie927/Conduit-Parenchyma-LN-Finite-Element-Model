function [Pressure_Sol] = Fluid_ransport_Solver(NL,EL,NoN,Ly,Line_list, width_x,width_y,Conduit_Surf, bound_node1, bound_node2,conduit_coordinates1,conduit_coordinates2,conduit_coordinates3,conduit_coordinates4,Line_list2)

%% Parameters % Use cgs unit system for all the parameters.
mu = 1E-02; % [(dyne s)/cm^2]  Lymph viscosity
k_p = 5E+4; % [cm^2] Permeability tensor
rou = 1000; % [kg/m^3] Fluid density
g = 9.81; % [m^2/s] Gravitational constant
K_cond = (mu/k_p)*0.1; % [] Pseudo hydraulic conductivity
f = 1000; % [cm^2/s] mass flow source term
dt = 1e-3; % Time-step [s]
total_time = 0.1; % Total time for PDE [s]
no_t_steps = total_time/dt; % Number of time-steps [ ]
% NoN = (nx+1)*(ny+1); % Total number of nodes 2D
L = 1e-2; % Characteristic Length
p9 = 1; p10 = 1; p11 = 1; p12 = 1; p15 = 1; p27 = 1; p99 = 0;
b1 = 0; b2 = 0; b3 = 0; b4 = 0; b5 = 0; b6_conduit = 0;
q = 0; % Flux on Left Boundary (Von Neumann condition)
Pb_left = 0; % Direchlet condition on Right Boundary
Pb_right = 0;
Pb_bottom = 0;
Pb_Middle = 0;
Pb_top = 0;
k = 1000000000; % Parameter in Direchlet BC % Flux of concentration.

%% Meshing 2D
% [xmesh ymesh] = meshgrid(linspace(1,width_x,n),linspace(1,width_y,n));
% [NL, EL] = meshing(width_x,width_y,nx,ny,'triangle');
% load("Conduit_CS_Mesh_segment_v2.mat")
% NL = node_Pos(:,1:2);
% EL = El_TRIANGLES_connectivity;
% NoN = size(node_Pos,1);
% Ly = width_y/(round(sqrt(NoN)));

%% Boundary Conditions
% Find Surfaces
[Left_Surf,Left_connectivity,Left_indices,Right_Surf,Right_connectivity,Right_indices,Bottom_Surf, Bottom_connectivity, Bottom_indices,Top_Surf, Top_connectivity,Top_indices] = FindSurfaces(NL,EL,width_x,width_y);

%% Solving PDE Weak Form Coefficients
% Find Parameter A, a, b, and c
% Find Vector for Area of Triangles
[A,a,b,c] = ShapeFunctionCoeff(NL,EL);

%% Pre-allocate concentration matrixes (and others)
% C = zeros(NoN,NoN); % Concentration values at node positions
array_size = size(EL,1);
diff_term = zeros(3,3);
source_term = zeros(3,1);
von_neum_cond_5th = zeros(3,1);
dir_term_3rd = zeros(3,3);
direchlet_cond_4th = zeros(3,1);
debug_1 = 0; debug_2 = 0; debug_3 = 0; debug_4 = 0; debug_5 = 0; peter = 0;

%% 1-D Lumped Parameter Model
R1 = 1E-5; % [Dynes/cm2/ul] % Poiseuille Resistance
R2 = 1E-5; % [Dynes/cm2/ul] % Poiseuille Resistance
R3 = 1E5;
k_p1 = (R2*R3)/(R1*R2+R1*R3+R2*R3);
k_p3 = (R1*R3)/(R1*R2+R1*R3+R2*R3);

C11 = (k_p1/R1)-1/R1;
C12 = k_p3/R1;
C13 = 0;
C21 = (k_p1/R2);
C22 = (k_p3/R2)-1/R2;
C23 = 0;
C31 = -k_p1/R3;
C32 = -k_p1/R3;
C33 = +1;

K_el1 = [C11 C12 C13;
        C21 C22 C23;
        C31 C32 C33]; % 1st node is boundary

K_el2 = [0 1/(2*R3) 1/(2*R3);
         +(R3*Ly)/(2*K_cond) -Ly/(4*K_cond) -Ly/(4*K_cond);
         +(R3*Ly)/(2*K_cond) -Ly/(4*K_cond) -Ly/(4*K_cond);];

params = [R1 R2 R3];
P_bc = [0 0]; P1 = P_bc(1); P2 = P_bc(2);

[U, Conduit_flow] = Fluid_LPM_assembly_v2(NL, EL, Conduit_Surf, params, bound_node1, bound_node2, P_bc);


Left_bound=sum(ismember(EL,Left_indices),2)==2;
Right_bound=sum(ismember(EL,Right_indices),2)==2;
Top_bound=sum(ismember(EL,Top_indices),2)==2;
Bottom_bound=sum(ismember(EL,Bottom_indices),2)==2;
Conduit_bound=sum(ismember(EL,Line_list2(1)),2)==1;
Leftidx = find(Left_bound); Rightidx = find(Right_bound);
Topidx = find(Top_bound); Bottomidx = find(Bottom_bound); 
Conduitidx = find(Conduit_bound);
[tf_left Leftidx]=ismember(EL(Leftidx,:), Left_indices);
[tf_right Rightidx]=ismember(EL(Rightidx,:), Right_indices);
[tf_top Topidx]=ismember(EL(Topidx,:), Top_indices);
[tf_bottom Bottomidx]=ismember(EL(Bottomidx,:), Bottom_indices);
[tf_conduit Conduitidx2]=ismember(EL(Conduitidx,:), Line_list2);
% map=bsxfun(@eq, EL(Leftidx),Left_indices(:));
bcon = unique(Conduit_Surf);
no_el = size(Conduit_Surf,1);
no_nodes = length(bcon); bcon2 = sort(bcon);
Conduit_Surf_pseudo = zeros(no_nodes,2);
NoLPM = size(Conduit_Surf,1) + no_nodes;

%% Initialize Matrices
K = zeros(NoN+NoLPM,NoN+NoLPM);
P = zeros(NoN+NoLPM,1);
%% 2D Galerkin Weak Form
for i=1:size(EL,1)  % solve Weak Form PDE for each element
    LE = flip(EL(i,1:2));
    if Left_bound(i)
        b1 = b1+1;
        loc_idx = find(tf_left(b1,:));
        % von_neum_cond_5th(loc_idx,1) = -q*(Ly/2)*[1;1]; %5th term
        dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        direchlet_cond_4th(loc_idx,1)=-((k*Pb_left*Ly)/2).*[1;1]; %4th term
        debug_1 = debug_1+1;
        Left_El(debug_1,:) = EL(i,:);
    elseif Right_bound(i) % why first two elements lie on edge?
        b2 = b2+1;
        loc_idx = find(tf_right(b2,:));
        dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        direchlet_cond_4th(loc_idx,1)=-((k*Pb_right*Ly)/2).*[1;1]; %4th term
        debug_2 = debug_2+1;
    elseif Top_bound(i)
        b3 = b3+1;
        loc_idx = find(tf_top(b3,:));
        dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        direchlet_cond_4th(loc_idx,1)=-((k*Pb_top*Ly)/2).*[1;1]; %4th term
        debug_3 = debug_3+1;
        % Top_El(debug_3,:) = EL(i,:);
    elseif Bottom_bound(i)
        b4 = b4+1;
        loc_idx = find(tf_bottom(b4,:));
        dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        direchlet_cond_4th(loc_idx,1)=-((k*Pb_bottom*Ly)/2).*[1;1]; %4th term
        debug_4 = debug_4+1;
    elseif Conduit_bound(i)
        % b5 = b5+1;
        % loc_idx = find(tf_conduit(b5,:));
        % dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        % direchlet_cond_4th(loc_idx,1)=-((k*P_bc(1)*Ly)/2).*[1;1]; %4th term
        % debug_4 = debug_4+1;
        %(ismember(EL(i,1),Line_list2(:,1))) && (ismember(EL(i,2),Line_list2(:,2)))  %(ismember(LE,Line_list2)); %&& (ismember(EL(i,2),Line_list(:,2))) %|| (ismember(EL(i,2),Line_list(conduit_indices(:),1)) && (ismember(EL(i,1),Line_list(conduit_indices(:),2))))
        % dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term %% Write a truth table for this depending on what node numbers are located in the element
        % find which edge the element belongs to
        % edge_no = find(ismember(EL(i,1:2),Line_list2,'rows'));
        % if isempty(edge_no) == 1
        %     edge_no = find(ismember(flip(EL(i,1:2)),Line_list2,'rows'));
        % direchlet_cond_4th(:,:)=-((k*(Conduit_flow(edge_no)*R3)*Ly)/2).*[1;1;0]; %4th term
        for j = 1:size(Line_list2,1)
            if (ismember(Line_list2(j,:),EL(i,:)))
                b6_conduit = j;
            end
        end
        % dir_term_3rd(loc_idx,loc_idx) = -k*(Ly/6)*[2,1;1,2];  %3rd term %% Write a truth table for this depending on what node numbers are located in the element
        % direchlet_cond_4th(loc_idx,1)=((k*(Conduit_flow(b6_conduit)*R3)*Ly)/2).*[1;1]; %4th term
        % von_neum_cond_5th(:,:) = ((Conduit_flow(b6_conduit)*R3))*(Ly/2)*[1;1;0]; %5th term
        p15 = p15+1;
    else
        diff_term(:,:) =(-1/(2*A(i))).*(b(i,:)'*b(i,:)+c(i,:)'*c(i,:));     %1st term
        % source_term(:,:) = 0;     %2nd term
    end
    %% Assemble matrices and vectors []
    for i2 = 1:3
        for i3 = 1:3
            nno = EL(i,:);      % Return the node numbers in the element
            dir_temp = dir_term_3rd(i2,i3);
            diff_temp = diff_term(i2,i3);
            K(NoLPM+nno(i2),NoLPM+nno(i3)) = K(NoLPM+nno(i2),NoLPM+nno(i3))+dir_temp+diff_temp;
        end
        vN_temp_5th = von_neum_cond_5th(i2,1);
        dir_temp_4th = direchlet_cond_4th(i2,1);
        adv_temp = source_term(i2,1);
        P(NoLPM+nno(i2),1)  = P(NoLPM+nno(i2),1) + dir_temp_4th + vN_temp_5th;
    end
    diff_term = zeros(3,3);
    source_term = zeros(3,1);
    dir_term_3rd = zeros(3,3);
    von_neum_cond_5th = zeros(3,1);
    direchlet_cond_4th = zeros(3,1);
end

for i = 1:no_el
    Conduit_Surf_pseudo(i,:) = [find(bcon==Conduit_Surf(i,1)) find(bcon==Conduit_Surf(i,2))];
end

for j = 1:length(bound_node1)
    [~,bound_node1_mapped(j)]=ismember(bound_node1(j),bcon);
end
for j = 1:length(bound_node2)
    [~,bound_node2_mapped(j)]=ismember(bound_node2(j),bcon);
end

for el_no = 1:no_el
    nno = Conduit_Surf_pseudo(el_no,:);
    nno_r = Conduit_Surf(el_no,:);
    nno_a_q = [nno no_nodes+el_no];
    nno_a_q2 = [no_nodes+el_no NoLPM+nno_r];
    K(nno_a_q,nno_a_q) = K(nno_a_q,nno_a_q) + K_el1;
    K(nno_a_q2,nno_a_q2) = K(nno_a_q2,nno_a_q2) + K_el2;
end

%% LPM boundary conditions
for i = 1:length(bound_node1) % Replace flow boundary condition with a pressure boundary condition
    flow_cons_v = zeros(1,NoN+NoLPM);
    P(bound_node1_mapped(i)) = P1;
    flow_cons_v(bound_node1_mapped(i)) = 1;
    K(bound_node1_mapped(i),:) = flow_cons_v;
end
for i = 1:length(bound_node2) % Replace flow boundary condition with a pressure boundary condition
    flow_cons_v = zeros(1,NoN+NoLPM);
    P(bound_node2_mapped(i)) = P2;
    flow_cons_v(bound_node2_mapped(i)) = 1;
    K(bound_node2_mapped(i),:) = flow_cons_v;
end

row_all_zeros = find(all(K == 0,2));
C = pinv(K)*P; % why are elements containing nodes of the edge defining mesh sometimes not included in the Element Edge List of the FEM?
%% reconstruct Q3 Flows
Flow_terms = C(((no_nodes+1):NoLPM),1);
Conduit_flow_v2 = zeros(no_el,1);
for el_no = 1:no_el
    nno_r = Conduit_Surf(el_no,:);
    Conduit_flow_v2(el_no) = Flow_terms(el_no) + (1/(2*R3))*(C(NoLPM+nno_r(1))+C(NoLPM+nno_r(2)));
end

Flow_val = (Conduit_flow_v2 - Conduit_flow);
Pressure_Sol = C((NoLPM+1):(NoLPM+NoN));
end
