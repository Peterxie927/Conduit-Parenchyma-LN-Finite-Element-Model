function [C4 C_ot dt total_time] = ADR_transport_solver(NL,EL,NoN,Ly,Pe,Line_list,width_x,width_y,Conduit_Surf, bound_node1, bound_node2,conduit_coordinates1,conduit_coordinates2,conduit_coordinates3,conduit_coordinates4,Line_list2)
%% Parameters
D = 1E-10; % Diffusion constant [m^2/s]
ux = 0E-14; % Velocity [m/s]
uy = 0E-14; % Velocity [m/s]
dt = 1E-05; % Time-step [s]
total_time = 1E-01; % Total time for PDE [s]
no_t_steps = total_time/dt; % Number of time-steps [ ]
% nx = 50; % Number of Nodes x-direction [ ]
% ny = 50; % Number of Nodes y-direction [ ]
n1D = 4;
NoN = size(NL,1); % Total number of nodes 2D
L = 1e-2; % Characteristic Length
p9 = 1; p10 = 1; p11 = 1; p12 = 1; p15 = 1; p27 = 1; p99 = 0;
b1 = 0; b2 = 0; b3 = 0; b4 = 0; b5 = 0; b6_conduit = 0;
Pe = [(ux*width_x/D) (uy*width_y/D)]; % Dimensionless Peclet Number. Oscillations in the Peclet number exist when the Peclet number exceeds the value of 10.
q = 0; % Flux on Left Boundary (Von Neumann condition)
Cb_left = 0; % Direchlet condition on Left Boundary
Cb_right = 0; % Direchlet condition on Right Boundary
Cb_bottom = 0; % Direchlet condition on Bottom Boundary
Cb_top = 0; % Direchlet condition on Top Boundary
k = 1000000; % Parameter in Direchlet BC % Flux of concentration.
debug_1 = 0; debug_2 = 0; debug_3 = 0; debug_4 = 0; debug_5 = 0; peter = 0; p_count = 0; p99 = 0;

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
diff_term = zeros(3,3,array_size);
adv_term = zeros(3,3,array_size);
von_neum_cond_5th = zeros(3,1,array_size);
dir_term_3rd = zeros(3,3,array_size);
direchlet_cond_4th = zeros(3,1,array_size);
debug_1 = 0; debug_2 = 0; debug_3 = 0; debug_4 = 0; debug_5 = 0;
K = zeros(NoN,NoN);
U = zeros(NoN,1);
m_el = (1/12)*A(1,1)*[2,1,1;1,2,1;1,1,2];
for i = 1:size(Right_indices,2)
    C_ot(Right_indices(i),1) = Cb_right;
end
for i = 1:size(Left_indices,2)
    C_ot(Left_indices(i),1) = Cb_left;
end

Left_bound=sum(ismember(EL,Left_indices),2)==2;
Right_bound=sum(ismember(EL,Right_indices),2)==2;
Top_bound=sum(ismember(EL,Top_indices),2)==2;
Bottom_bound=sum(ismember(EL,Bottom_indices),2)==2;
Conduit_bound=sum(ismember(EL,Line_list2),2)==2;
Leftidx = find(Left_bound); Rightidx = find(Right_bound);
Topidx = find(Top_bound); Bottomidx = find(Bottom_bound);
Conduitidx = find(Conduit_bound);
[tf_left Leftidx]=ismember(EL(Leftidx,:), Left_indices);
[tf_right Rightidx]=ismember(EL(Rightidx,:), Right_indices);
[tf_top Topidx]=ismember(EL(Topidx,:), Top_indices);
[tf_bottom Bottomidx]=ismember(EL(Bottomidx,:), Bottom_indices);
[tf_conduit Conduitidx2]=ismember(EL(Conduitidx,:), Line_list2);
bcon = unique(Conduit_Surf);
no_el = size(Conduit_Surf,1);
no_nodes = length(bcon); bcon2 = sort(bcon);

for i=1:size(EL,1)  % solve Weak Form PDE for each element
    if Left_bound(i)
        % von_neum_cond_5th(loc_idx,loc_idx,i) = -q*(Ly/2)*[1;1]; %5th term
        b1 = b1+1;
        loc_idx = find(tf_left(b1,:));
        dir_term_3rd(loc_idx,loc_idx,i) = -k*(Ly/6)*[2,1;1,2];  %3rd term
        direchlet_cond_4th(loc_idx,1,i)=-((k*Cb_left*Ly)/2).*[1;1]; %4th term
        debug_1 = debug_1+1;
        Left_El(debug_1,:) = EL(i,:);
    elseif Right_bound(i)
        b2 = b2+1;
        loc_idx = find(tf_right(b2,:));
        dir_term_3rd(loc_idx,loc_idx,i) = -k*(Ly/6)*[2,1;1,2];   %3rd term
        direchlet_cond_4th(loc_idx,1,i)= -((k*Cb_right*Ly)/2).*[1;1]; %4th term
        debug_2 = debug_2+1;
    elseif Top_bound(i)
        b3 = b3+1;
        loc_idx = find(tf_top(b3,:));
        dir_term_3rd(loc_idx,loc_idx,i) = -k*(Ly/6)*[2,1;1,2];   %3rd term
        direchlet_cond_4th(loc_idx,1,i)=-((k*Cb_top*Ly)/2).*[1;1]; %4th term
        debug_3 = debug_3+1;
        Top_El(debug_3,:) = EL(i,:);
    elseif Bottom_bound(i)
        b4 = b4+1;
        loc_idx = find(tf_bottom(b4,:));
        dir_term_3rd(loc_idx,loc_idx,i) = -k*(Ly/6)*[2,1;1,2]; %3rd term
        direchlet_cond_4th(loc_idx,1,i)=-((k*Cb_bottom*Ly)/2).*[1;1]; %4th term
        debug_4 = debug_4+1;
        Bottom_El(debug_4,:) = EL(i,:);
        % elseif Conduit_bound(i)
        %     b5 = b5+1;
        %     loc_idx = find(tf_conduit(b5,:));
        %     for j = 1:size(Line_list2,1)
        %         if (ismember(Line_list2(j,:),EL(i,:)))
        %             b6_conduit = j;
        %         end
        %     end
    else
        diff_term(:,:,i) =-(1/(2*A(i))).*(b(i,:)'*b(i,:)+c(i,:)'*c(i,:)); %1st term
        adv_term(:,:,i) = -(1/3)*Pe(1)*[b(i,:);b(i,:);b(i,:)]+Pe(2)*[c(i,:);c(i,:);c(i,:)]; %2nd term
        if Conduit_bound(i)
            b5 = b5+1;
            loc_idx = find(tf_conduit(b5,:));
            for j = 1:size(Conduit_Surf,1)
                if (ismember(Conduit_Surf(j,:),EL(i,:)))
                    b6_conduit = j;
                    % insert flux terms here
                end
            end

        end
    end
end

%% Place the 2D FDM Model as a function in here
%% Finite Difference Function
Dr = 1E-14; Dz = 2E-10; alpha = 1; ur = 0E-05; uz = 0e-015;
[C4, K_FDM, out_flux, border_C, dr, left_index, right_index, NL_pseudo1,NL_pseudo2,Conduit_Surf_pseudo] = ADR_FD_Network_function(NL,EL,Conduit_Surf, bound_node1, bound_node2, Ly, Dr, Dz, ur, uz, alpha,dt);

Nr = 6; No2DFDM = no_nodes*(Nr+1);
EL_connectivity = Conduit_Surf_pseudo;
EL_connectivity2 = [];
for p = 1:size(EL_connectivity,1)
    EL_connectivity2 = [EL_connectivity2 EL_connectivity(p,:)];
end
bcon = unique(EL_connectivity);
no_nodes = length(bcon);
node_counts = histcounts(EL_connectivity);  % algorithm to enforce boundary condition for nodes that appear once only (assume all of them are sinks with C=0)
boundary_nodes = find(node_counts == 1);
z_coords = zeros(1,no_nodes);
Nz = size(Conduit_Surf_pseudo,1)-1;
%% Boundary Conditions
Ct=0;                              % r = 0 Dirichlet B.C
Cb=0;                              % r = Lr Dirichlet B.C
Cl=0;                              % z = 0 Dirichlet B.C (Source Nodes concentration)
Cr=0;                              % z = Lz Dirichlet B.C (Boundary Nodes concentration)
CvNt=0*dr;                              % r = 0 von Neumann B.C IMPEMENT THIS ZERO-FLUX BOUNDARY CONDITION
d_mem = 1E-06;
% CvNb=0;                              % r = Lr von Neumann B.C
% CvNl=0;                              % z = 0 von Neumann B.C
% CvNr=0;                              % z = Lz von Neumann B.C
%% Initialize Matrices
K2 = zeros(NoN+No2DFDM,NoN+No2DFDM);
U2 = zeros(NoN+No2DFDM,1);

K = zeros(NoN+No2DFDM,NoN+No2DFDM);
U = zeros(NoN+No2DFDM,1);
M = zeros(NoN+No2DFDM,NoN+No2DFDM);
C = zeros(NoN+No2DFDM,1);
C_ot = zeros(NoN+No2DFDM,no_t_steps);

%% Assemble matrices and vectors []
for i = 1:size(EL,1)
    for i2 = 1:size(EL,2)
        for i3 = 1:size(EL,2)
            nno = EL(i,:);      % Return the node numbers in the element
            dir_temp_3rd = dir_term_3rd(i2,i3,i);
            diff_temp = diff_term(i2,i3,i);
            adv_temp = adv_term(i2,i3,i);
            K(nno(i2),nno(i3)) = K(nno(i2),nno(i3))+dir_temp_3rd+diff_temp+adv_temp;
            M(nno(i2),nno(i3)) = M(nno(i2),nno(i3))+m_el(i2,i3);
        end
        vN_temp_5th = von_neum_cond_5th(i2,1,i);
        dir_temp_4th = direchlet_cond_4th(i2,1,i);
        U(nno(i2),1)  = U(nno(i2),1) + dir_temp_4th + vN_temp_5th;
    end
end
K((NoN+1):(NoN+No2DFDM),(NoN+1):(NoN+No2DFDM)) = -eye(No2DFDM);
M((NoN+1):(NoN+No2DFDM),(NoN+1):(NoN+No2DFDM)) = eye(No2DFDM);
for j = 1:length(bound_node1)
    [~,bound_node1_mapped(j)]=ismember(bound_node1(j),b);
end
for j = 1:length(bound_node2)
    [~,bound_node2_mapped(j)]=ismember(bound_node2(j),b);
end

C_fd = (alpha*Ly)/(4*dr);
C_fe = -(alpha*Ly)/(4*dr);

K_el1 = [0 0 0 0
    0 0 0 0
    C_fd C_fd C_fe C_fe
    C_fd C_fd C_fe C_fe];
% K_el1 = [0 0 0 0
%     0 0 0 0
%     C_fd 0 C_fe 0
%     C_fd 0 C_fe 0];
% K_el1 = [C_fd;
%          C_fd;
%          C_fe;
%          C_fe];
% for el_no = 1:no_el
%     nno = Conduit_Surf_pseudo(el_no,:);
%     nno_r = Conduit_Surf(el_no,:);
%     nno_a_q = [NoN+nno nno_r];
%     K(nno_a_q,nno_a_q) = K(nno_a_q,nno_a_q) + K_el1;
%     K(nno_a_q2,nno_a_q2) = K(nno_a_q2,nno_a_q2) + K_el2;
% end

%%
%   % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 2;
%   % set the seconds per image
% % open the video writer
%   open(writerObj);

t_real = 0;
C_ot((NoN+1):(NoN+No2DFDM),1) = C4;
peter = figure('Renderer', 'painters', 'Position', [10 10 1000 700])

    for el_no = 1:no_el
        nno = Conduit_Surf_pseudo(el_no,:);
        nno_r = Conduit_Surf(el_no,:);
        nno_a_q = [NoN+(nno)*(Nr+1) nno_r];
        K(nno_a_q,nno_a_q) = K(nno_a_q,nno_a_q) + K_el1;
    end
    
for t = 2:no_t_steps
    t_real = t_real+dt;
    % C_ot((NoN+1):(NoN+No2DFDM),t) = C4;
    C5 = C4;
    vN_c = zeros(no_nodes,1);
   
    Peter = pinv(M+dt*(-K));
    C_ot((NoN+1):(NoN+No2DFDM),t) = C4;
    C_ot(:,t) = Peter*(M*C_ot(:,(t-1))+dt*(-U));    
    C_ot((NoN+1):(NoN+No2DFDM),t) = C4;

    for node_no = 1:no_nodes
        % nno_r = Conduit_Surf(el_no,:);
        % nno = Conduit_Surf_pseudo(el_no,:);
        vN_c(node_no,1) = alpha*(C_ot(bcon2(node_no),t-1)-C4(bcon(node_no)*(Nr+1)));
    end

    C5(left_index)=C5(left_index+1)-CvNt*dr; % Zero flux boundary condition at r=0;
    C5(right_index)=C5(right_index-1)+vN_c*dr;  % Zero flux boundary condition at r=0;
    C4 = K_FDM*C5;

    if t == 1000
        P_title = sprintf('t = %f',t_real+dt);
        title(P_title);
        subplot(2,2,1)
        patch('Faces',EL,'Vertices',NL,'CData',C_ot(1:NoN,t),'FaceColor','interp');
        axis equal;
        box on; grid on;
        axis([0 width_x 0 width_y])
        hold on
        for i = 1:size(conduit_coordinates1,1)
            plot([conduit_coordinates1(i) conduit_coordinates3(i)],[conduit_coordinates2(i) conduit_coordinates4(i)],'LineWidth',1.8,'color','k')
            hold on
        end
        set(gcf,'color','w'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
        hold off
        % P_title = sprintf('t = %f',t_real);
        % title(P_title);
        % xlabel('X [ ]','FontSize',12,'interpreter','Latex');
        % ylabel('Y [ ]','FontSize',12,'interpreter','Latex');
        colorbar;
        % F(i) = getframe(gcf);
        subplot(2,2,2)
        for p = 1:size(EL_connectivity,1)
            patch([NL_pseudo1(p,1) NL_pseudo2(p,1)],[NL_pseudo1(p,2) NL_pseudo2(p,2)],[z_coords(EL_connectivity(p,:))],[C4(right_index(EL_connectivity(p,:)))],'facecolor','none','edgecolor','interp','linewidth',2);
            hold on
        end
        hold off
        % shading interp
        colorbar;colormap(jet);
        % xlabel('X [ ]','FontSize',16,'interpreter','Latex');
        % ylabel('Y [ ]','FontSize',16,'interpreter','Latex');
        axis([0 width_x 0 width_y]);
        p_count = p_count+1;
        % peter_title = sprintf('2D ADR Dr = %d, Dz = %d, ur = %d, uz = %d', Dr, Dz, ur, uz);
        % peter_title = sprintf('Coupling alpha = %d',alpha);
        % title({peter_title;['time = ',num2str(t_real)]})
        set(gcf,'color','w'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        box on; grid on;
        drawnow;
    end
     % for i = 1:length(bound_node1_mapped)
     %        C4(left_index(bound_node1_mapped(i)):(left_index(bound_node1_mapped(i))+Nr)) = Cl;
     % end
        % for i = 1:length(bound_node2_mapped)
        %     C4(left_index(bound_node2_mapped(i)):(left_index(bound_node2_mapped(i))+Nr)) = Cr;
        % end
        % % C4(right_index) = Cb;
        % C4(left_index) = Ct;
        % % Neumann BCs:
        % C4(left_index)=C4(left_index+1)-CvNt*dr; % Zero flux boundary condition at r=0;
    % end
    p99 = p99+1;
    % F(p99) = getframe(gcf);
    disp(t);
    toc
end

% writerObj = VideoWriter('Conduit_Network_Coupling_v3_alpha1.avi');
% writerObj.FrameRate = 100;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:(length(F)-402)
%     % convert the image to a frame
%     frame = F(i) ;
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

end

