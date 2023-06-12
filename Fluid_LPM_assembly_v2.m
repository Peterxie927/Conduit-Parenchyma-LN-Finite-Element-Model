function [U, Conduit_flow] = Fluid_LPM_assembly(NL,EL,Conduit_Surf, params, bound_node1, bound_node2, P_bc);

R1 = params(1); R2 = params(2); R3 = params(3);
k_p1 = (R2*R3)/(R1*R2+R1*R3+R2*R3);
k_p3 = (R1*R3)/(R1*R2+R1*R3+R2*R3);
i2 = 0; i3 = 0; i4 = 0;

C11 = (k_p1/R1)-1/R1;
C12 = k_p3/R1;
C13 = 0;
C21 = (k_p1/R2);
C22 = (k_p3/R2)-1/R2;
C23 = 0;
C31 = -k_p1/R3;
C32 = -k_p1/R3;
C33 = -1;

K_el1 = [C11 C12 C13;
        C21 C22 C23;
        C31 C32 C33]; % 1st node is boundary



n_el = size(Conduit_Surf,1); % size of the 0D LPM 
b = unique(Conduit_Surf);
no_nodes = length(b); b2 = sort(b);
Conduit_Surf_pseudo = zeros(no_nodes,2);
kglobal = zeros(no_nodes+n_el,no_nodes+n_el);
U = zeros(no_nodes+n_el, 1);
f_lpm = zeros(no_nodes+n_el,1);

for i = 1:n_el
    Conduit_Surf_pseudo(i,:) = [find(b==Conduit_Surf(i,1)) find(b==Conduit_Surf(i,2))];
end

for i = 1:no_nodes
    if (sum(Conduit_Surf(:)==i) == 1) % boundary node
        i2 = i2+1;
        bound_node(i2,1) = i;
    else
        i3 = i3+1;
        unbound_node(i3,1) = i;
    end
end

for j = 1:length(bound_node1)
    [~,bound_node1_mapped(j)]=ismember(bound_node1(j),b);
end
for j = 1:length(bound_node2)
    [~,bound_node2_mapped(j)]=ismember(bound_node2(j),b);
end

no_bound = length(bound_node1)+length(bound_node2);
no_unbound = length(unbound_node);
no_unknown = (no_nodes+n_el)-no_bound;

for el_no = 1:n_el
    nno = Conduit_Surf_pseudo(el_no,:);
    flow_cons_v = zeros(1,no_nodes+n_el);
    nno_a_q = [nno no_nodes + el_no];
        kglobal(nno_a_q,nno_a_q) = kglobal(nno_a_q,nno_a_q) + K_el1;
    if (ismember(nno(1),bound_node)) || (ismember(nno(2),bound_node)) 
        i4 = i4+1; % Replace flow boundary condition with a pressure boundary condition
    end
end

%% Flow Equations
rem_eq = no_unknown - el_no;

%% Boundary Conditions

P1 = P_bc(1); P2 = P_bc(2);
f_lpm = zeros(no_nodes+n_el,1);
Conduit_flow = zeros(n_el,1);
no_bound = length(bound_node1) + length(bound_node2);
% Map Pseudo Node numbering into actual node numbering

for i = 1:length(bound_node1) % Replace flow boundary condition with a pressure boundary condition
    flow_cons_v = zeros(1,no_nodes+n_el);
    f_lpm(bound_node1_mapped(i)) = P1;
    flow_cons_v(bound_node1_mapped(i)) = 1;
    kglobal(bound_node1_mapped(i),:) = flow_cons_v;
end
for i = 1:length(bound_node2) % Replace flow boundary condition with a pressure boundary condition
    flow_cons_v = zeros(1,no_nodes+n_el);
    f_lpm(bound_node2_mapped(i)) = P2;
    flow_cons_v(bound_node2_mapped(i)) = 1;
    kglobal(bound_node2_mapped(i),:) = flow_cons_v;
end

U = kglobal\f_lpm;
Conduit_flow = U(((no_nodes+n_el)-n_el+1):(no_nodes+n_el));
end

