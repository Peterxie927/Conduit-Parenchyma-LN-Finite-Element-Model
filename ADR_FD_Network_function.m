function [C4, K, out_flux, border_C, dr, left_index,right_index,NL_pseudo1,NL_pseudo2,Conduit_Surf_pseudo] = ADR_FD_Network_function(NL,EL,Conduit_Surf, bound_node1, bound_node2, Ly, Dr, Dz, ur, uz, alpha,dt)

NL = NL(:,1:2);
EL = EL(:,1:3);
NoN = size(NL,1);
width_y = 1; width_x = 1; p27 = 1; p11 = 1; p10 = 1; i2 = 0; i3 = 0;
Ly = width_y/(round(sqrt(NoN)));

n_el = size(Conduit_Surf,1); % size of the 0D LPM
b = unique(Conduit_Surf);
no_nodes = length(b); b2 = sort(b);
Conduit_Surf_pseudo = zeros(n_el,2);

for i = 1:n_el
    Conduit_Surf_pseudo(i,:) = [find(b==Conduit_Surf(i,1)) find(b==Conduit_Surf(i,2))];
    NL_pseudo1(i,:) = [NL(b(Conduit_Surf_pseudo(i,1)),1) NL(b(Conduit_Surf_pseudo(i,1)),2)];
    NL_pseudo2(i,:) = [NL(b(Conduit_Surf_pseudo(i,2)),1) NL(b(Conduit_Surf_pseudo(i,2)),2)];
end

EL_connectivity = Conduit_Surf_pseudo;
EL_connectivity2 = [];
for p = 1:size(EL_connectivity,1)
    EL_connectivity2 = [EL_connectivity2 EL_connectivity(p,:)];
end
bcon = unique(EL_connectivity);
no_nodes = length(bcon);
node_counts = histcounts(EL_connectivity);  % algorithm to enforce boundary condition for nodes that appear once only (assume all of them are sinks with C=0)
boundary_nodes = find(node_counts == 1);
% x_coords = [1.45859683191491	2.14291950941019	3.37860011455536	4.37686454713925	5.19022292348768	6.28391082036261	7.03792714478153	8.02697505933330]*1E-06;
% y_coords = [1.26539877650449	1.38958361505101	1.46700534211459	1.06495310423687	1.28441183043610	1.23469532052910	1.00595103475062	1.16856132219944]*1E-06;

% for i = 1:1:size(bcon,1)
%     x_coords(i) = NL(b(i),1);
%     y_coords(i) = NL(b(i),2);
% end
z_coords = zeros(1,no_nodes);
[X, Y] = meshgrid(1:no_nodes, 1:no_nodes);

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


Lr = 1E-05;                      % [m] Conduit radius
Lz = 1E-05;                      % [m] Conduit segment length
% EL_connectivity = [1 2; 2 3; 3 4;3 8; 4 5; 5 6;5 7;]; % 4 segments
% EL_connectivity = [1 2; 2 3; 3 4;4 5;5 6;6 7;7 8;];


% Add a random perturbation to the x and y coordinates
% x_coords = X(1,:) + rand(1,no_nodes)*0.5;
% y_coords = Y(1,:) + rand(1,no_nodes)*0.5;
% x_coords = linspace(1E-06,1E-04,no_nodes)+rand(1,no_nodes)*5E-06;
% y_coords = linspace(1E-06,1E-04,no_nodes)+rand(1,no_nodes)*5E-06;
[X, Y] = meshgrid(linspace(0,1,100), linspace(0,1,100));
Nz = no_nodes-1;                 % [ ] Number of z segments. Nodes in z direction = Nz+1
Nr = 6;                          % [ ] Number of r segments. Nodes in r direction = Nr+1
Dr = Dr;                      % [m^2/s] Diffusion Coefficient radially
Dz = Dz;                      % [m^2/s] Diffusion Coefficient axially
% t_end = 1E-3;                  % [s] Total simulation time
% dt = 1E-08;                    % [s] Time step (1 sec)
alpha = alpha;
dt = dt; % dt = 1E-06;                     % [s] Time step (1 sec)
dr = Lr/(Nr);                    % [m] delta r
dz = Lz/Nz;                    % [m] delta z
x = 0:dr:Lr;                     % Range of x(0,Lz) and specifying the grid points
y = 0:dz:Lz;                     % Range of y(0,Lr) and specifying the grid points
xp = x.*linspace(1,size(x,2),size(x,2));
mass_balance_vec = zeros();
ur = ur; uz = uz;
Pe_r_local = (ur*dr)/Dr; Pe_z_local = (uz*dz)/Dz;
Pe_r_global = (ur*Lr)/Dr; Pe_z_global = (uz*Lr)/Dz; % Peclet numbers between conduits and surrounding parenchyma?
stability_t_thres = (dz^2)/(uz*dz+2*Dz);
NoN = (Nr+1)*(Nz+1);
C4 = zeros(NoN,1); C5 = zeros(NoN,1);
K = zeros(NoN,NoN);
K2 = zeros(NoN,NoN); K3 = zeros(NoN,NoN);
p_count = 0;

%% Boundary Conditions
Ct=0;                              % r = 0 Dirichlet B.C
Cb=0;                              % r = Lr Dirichlet B.C
Cl=1;                              % z = 0 Dirichlet B.C (Source Nodes concentration)
Cr=0;                              % z = Lz Dirichlet B.C (Boundary Nodes concentration)
for i10 = 0:(Nz)
    C_nodes(i10+1) = 1+i10*(Nr+1);
end
CvNt=0*dr;                              % r = 0 von Neumann B.C IMPEMENT THIS ZERO-FLUX BOUNDARY CONDITION
% CvNb=0;                              % r = Lr von Neumann B.C
% CvNl=0;                              % z = 0 von Neumann B.C
% CvNr=0;                              % z = Lz von Neumann B.C

%% Unit Matrices Implementation and Boundary conditions
tic
left_index = 1:(Nr+1):NoN; right_index = (Nr+1):(Nr+1):NoN; top_index = 1:1:(Nr+1); bottom_index = (NoN-Nr):1:NoN;
t1 = (Dr*(dt/dr^2))-(Dr*(dt./(2.*x(2:(Nr+1)).*dr)))+((ur*dt)/(2*dr));
% t1_alpha = alpha*((Dr*(dt/dr^2))-(Dr*(dt/(2*x((Nr+1))*dr))))+((ur*dt)/(2*dr));
t2 = (Dr*(dt/dr^2))+(Dr*(dt./(2.*x(2:(Nr+1)).*dr)))-((ur*dt)/(2*dr));
t3 = 1-(2*Dr*dt)/(dr^2)-(2*Dz*dt)/(dz^2);
% t33 = 1-(2*Dr*dt)/(dr^2)-(1*Dz*dt)/(dz^2);
t32 = 1-(2*Dz*dt)/(dz^2);
t4 = (Dz*dt)/(dz^2) + (uz*dt)/(2*dz);
t5 = (Dz*dt)/(dz^2) - (uz*dt)/(2*dz);

% t1 = ((dt/dr^2))-((dt./(2.*x(2:(Nr+1)).*dr)))+((Pe_r_global*dt)/(2*dr));
% t2 = ((dt/dr^2))+((dt./(2.*x(2:(Nr+1)).*dr)))-((Pe_r_global*dt)/(2*dr));
% t3 = 1-(2*dt)/(dr^2)-(2*dt)/(dz^2);
% % t33 = 1-(2*Dr*dt)/(dr^2)-(1*Dz*dt)/(dz^2);
% t32 = 1-(2*dt)/(dz^2);
% t4 = (dt)/(dz^2) + (Pe_r_global*dt)/(2*dz);
% t5 = (dt)/(dz^2) - (Pe_z_global*dt)/(2*dz);

%% inlet BC
for i = 1:size(bound_node1_mapped,2)
    C4(left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr)) = Cl;
end
for i = 1:size(bound_node2_mapped,2)
    C4(left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr)) = Cr;
end
% C4(right_index) = Cb;
C4(left_index) = Ct;
% C4(1) = max(Ct,Cl); C4(Nr+1)=max(Cb,Cl);
% C4(NoN)=max(Cb,Cr);
% C4(NoN-Nr)=max(Ct,Cr);
% Neumann BCs:
C4(left_index)=C4(left_index+1)-CvNt*dr; % Zero flux boundary condition at r=0;

cent_mat = zeros(Nr+1,Nr+1);
prev_mat = zeros(Nr+1,Nr+1);
next_mat = zeros(Nr+1,Nr+1);
bc_mat = eye(Nr+1);
co1 = 0; co2 = 0;

%% Assemble Global stiffness matrix
for i = 1:(Nr+1)
    cent_mat(i,i) = t3;
    prev_mat(i,i) = t4;
    next_mat(i,i) = t5;
    if (i ~= (Nr+1)) && (i ~= (1))
        co1 = co1+1;
        cent_mat(i,i+1) = t2(co1);
    end
    if i ~= 1
        co2 = co2+1;
        cent_mat(i,i-1) = t1(co2);
    end
end

i=2:(Nr);
j=2:(Nz); % This index restricts allowing concentration at the edges to be free

[X,Y] = meshgrid(x,y);
P = [X(:) Y(:)];

%% Indices for matrices
x_nextflux = unique(EL_connectivity(:,1));
x_prevflux = unique(EL_connectivity(:,2));
N_nextflux = numel(x_nextflux);
N_prevflux = numel(x_prevflux);
count_next = zeros(N_nextflux,1);
count_prev = zeros(N_prevflux,1);
for k = 1:N_nextflux
    count_next(k) = sum(EL_connectivity(:,1)==x_nextflux(k));
end

for k = 1:N_prevflux
    count_prev(k) = sum(EL_connectivity(:,2)==x_prevflux(k));
end


    % 1st node: (1:Nr,1:Nr). 2nd node: ((Nr+1):2*Nr,(Nr+1):2*Nr). 3rd node:
    for co3 = 1:(no_nodes)
        K(((co3-1)*(Nr+1)+1):(co3*(Nr+1)),((co3-1)*(Nr+1)+1):(co3*(Nr+1))) = cent_mat;
        K(((co3-1)*(Nr+1)+1),((co3-1)*(Nr+1)+1):(co3*(Nr+1))) = 0;  K(((co3-1)*(Nr+1)+1),((co3-1)*(Nr+1)+1)) = t32; % only z-dependent for the first 'c' term
    end
    for co4 = 1:(size(EL_connectivity,1))
        % K((((EL_connectivity(co4,1)-1)*(Nr+1))+1):((EL_connectivity(co4,1))*(Nr+1)),(((EL_connectivity(co4,2)-1)*(Nr+1))+1):((EL_connectivity(co4,2))*(Nr+1))) = next_mat;
        % K((((EL_connectivity(co4,2)-1)*(Nr+1))+1):((EL_con\nectivity(co4,2))*(Nr+1)),(((EL_connectivity(co4,1)-1)*(Nr+1))+1):((EL_connectivity(co4,1))*(Nr+1))) = prev_mat;
        flux_divider = sum(EL_connectivity(:,1)==EL_connectivity(co4,1));
        K(left_index(EL_connectivity(co4,1)):left_index(EL_connectivity(co4,1))+Nr, left_index(EL_connectivity(co4,2)):left_index(EL_connectivity(co4,2))+Nr) = next_mat/flux_divider;
        K(left_index(EL_connectivity(co4,2)):left_index(EL_connectivity(co4,2))+Nr, left_index(EL_connectivity(co4,1)):left_index(EL_connectivity(co4,1))+Nr) = prev_mat;
    end

    for i = 1:length(bound_node1_mapped)
        K(left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr),1:NoN) = 0;
        K(left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr),left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr)) = bc_mat;
    end
    for i = 1:length(bound_node2_mapped)
        K(left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr),1:NoN) = 0;
        K(left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr),left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr)) = bc_mat;
    end
    K(1:(Nr+1),1:(Nr+1)) = bc_mat;

    for i = 1:length(bound_node1_mapped)
        C4(left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr)) = Cl;
    end
    for i = 1:length(bound_node2_mapped)
        C4(left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr)) = Cr;
    end
    % C4(right_index) = Cb;
    C4(left_index) = Ct;
    % C4(1) = max(Ct,Cl); C4(Nr+1)=max(Cb,Cl);
    % C4(NoN)=max(Cb,Cr);
    % C4(NoN-Nr)=max(Ct,Cr);
    % Neumann BCs:
    C4(left_index)=C4(left_index+1)-CvNt*dr; % Zero flux boundary condition at r=0;
    % C4(nx,:)=C4(nx-1,:)+CvNb*dr;
    % C4(:,1)=C4(:,2)-CvNr*dz;
    % C4(:,ny)=C4(:,ny-1)+CvNl*dz;
    C_M = reshape(C4,Nr+1,Nz+1);
    out_flux = (C4(right_index)-C4(right_index-1))./dr;
    border_C = C4(right_index);

    % mass_diff_v = mass_balance(C_M,Dz,dz,Dr,dr,uz,ur);
    % mass_diff_t(it) = sum(mass_diff_v,'all')/numel(mass_diff_v);

% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 10;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

end

