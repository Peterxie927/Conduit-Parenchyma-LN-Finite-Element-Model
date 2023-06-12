function [Conduit_Surf C_ot] = ADR_FEM_Solver(NL,EL,NoN,Ly,Line_list)
%% Parameters
D = 1E-10; % Diffusion constant [m^2/s]
ux = 0E-14; % Velocity [m/s]
uy = 0E-14; % Velocity [m/s]
dt = 0.001; % Time-step [s]
total_time = 0.1; % Total time for PDE [s]
no_t_steps = total_time/dt; % Number of time-steps [ ]
% nx = 50; % Number of Nodes x-direction [ ]
% ny = 50; % Number of Nodes y-direction [ ]
n1D = 4;
width_x = 1; % Width of domain (in x-direction) [m]
width_y = 1; % Width of domain (in y-direction) [m]
NoN = size(NL,1); % Total number of nodes 2D
L = 1e-2; % Characteristic Length
p9 = 1; p10 = 1; p11 = 1; p12 = 1; p15 = 1; p27 = 1; p99 = 0;
b1 = 0; b2 = 0; b3 = 0; b4 = 0; b5 = 0; b6_conduit = 0;
Pe = [(ux*width_x/D) (uy*width_y/D)]; % Dimensionless Peclet Number. Oscillations in the Peclet number exist when the Peclet number exceeds the value of 10.
q = 0; % Flux on Left Boundary (Von Neumann condition)
Cb_left = 0; % Direchlet condition on Left Boundary
Cb_right = 100; % Direchlet condition on Right Boundary
Cb_bottom = 0; % Direchlet condition on Bottom Boundary
Cb_top = 0; % Direchlet condition on Top Boundary
k = 100000; % Parameter in Direchlet BC % Flux of concentration.
debug_1 = 0; debug_2 = 0; debug_3 = 0; debug_4 = 0; debug_5 = 0; peter = 0;

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
m_el = (1/12)*A(1,1)*[2,1,1;1,2,1;1,1,2]; %% Only for equilateral triangles. Need to calculate for individual elements
M = zeros(NoN,NoN);
C = zeros(NoN,1);
C_ot = zeros(NoN,no_t_steps);
for i = 1:size(Right_indices,2)
    C_ot(Right_indices(i),1) = Cb_right;
end
for i = 1:size(Left_indices,2)
    C_ot(Left_indices(i),1) = Cb_left;
end

figure (3)
h=colorbar;
xlabel('x [m]','FontSize',16,'interpreter','Latex');
ylabel('y [m]','FontSize',16,'interpreter','Latex');
set(gcf,'color','w')

%% Find nodes corresponding to the Conduit Network
for p9 = 1:size(Line_list,1)
    if ((NL(Line_list(p9,1),1)==0) && (NL(Line_list(p9,2),1) ==0)) || ((NL(Line_list(p9,1),1)==width_x) && (NL(Line_list(p9,2),1) ==width_x)) || ((NL(Line_list(p9,1),2)==0) && (NL(Line_list(p9,2),2) ==0)) || ((NL(Line_list(p9,1),2)==width_y) && (NL(Line_list(p9,2),2) == width_y))
        % ((ismember(Line_list(1),Left_indices)&&ismember(Line_list(2),Left_indices))||(ismember(Line_list(1),Right_indices)&&ismember(Line_list(2),Right_indices))||(ismember(Line_list(1),Top_indices)&&ismember(Line_list(2),Top_indices))||(ismember(Line_list(1),Bottom_indices)&&ismember(Line_list(2),Bottom_indices)))
        peter = peter+1;
    else
        Conduit_Surf(p27,:) = [Line_list(p9,1) Line_list(p9,2)];
        conduit_coordinates1(p27,:) = [NL(Line_list(p9,1),1)];
        conduit_coordinates2(p27,:) = [NL(Line_list(p9,1),2)];
        conduit_coordinates3(p27,:) = [NL(Line_list(p9,2),1)];
        conduit_coordinates4(p27,:) = [NL(Line_list(p9,2),2)];
        conduit_indices(p27) = p9;
        Line_list2(p27,:) = [Line_list(p9,1) Line_list(p9,2)];
        % bound_node1 are given source BCs, bound_node2 are given sink BCs.
        if (((NL(Line_list(p9,1),1)==0) && ((NL(Line_list(p9,1),2)==0.2))))
            bound_node1(p10) = [Line_list(p9,1)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==0.2)))
            bound_node1(p10) = [Line_list(p9,2)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,1),1)==0.2) && ((NL(Line_list(p9,1),2)==0)))
            bound_node1(p10) = [Line_list(p9,1)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,2),1)==0.2) && ((NL(Line_list(p9,2),2)==0)))
            bound_node1(p10) = [Line_list(p9,2)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,1),1)==0) && ((NL(Line_list(p9,1),2)==0.8)))
            bound_node1(p10) = [Line_list(p9,1)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==0.8)))
            bound_node1(p10) = [Line_list(p9,2)];
            p10 = p10+1;
        elseif ((NL(Line_list(p9,1),1)==0.2) && ((NL(Line_list(p9,1),2)==1)))
            bound_node1(p10) = [Line_list(p9,1)];
            p10=p10+1;
        elseif ((NL(Line_list(p9,2),1)==0.2) && ((NL(Line_list(p9,2),2)==1)))
            bound_node1(p10) = [Line_list(p9,2)];
            p10=p10+1;
        elseif ((NL(Line_list(p9,1),1)==0) && ((NL(Line_list(p9,1),2)==0.5))) %|| ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==0.5)))
            bound_node1(p10) = [Line_list(p9,1)];
            p10=p10+1;
            % elseif ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==0.5))) %|| ((NL(Line_list(p9,2),1)==0) && ((NL(Line_list(p9,2),2)==0.5)))
            %     bound_node2(p11) = [Line_list(p9,2)];
            %     p11=p11+1;
        elseif ((NL(Line_list(p9,1),1)==0.8) && ((NL(Line_list(p9,1),2)==1)))
            bound_node2(p11) = [Line_list(p9,1)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,2),1)==0.8) && ((NL(Line_list(p9,2),2)==1)))
            bound_node2(p11) = [Line_list(p9,2)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,1),1)==1) && ((NL(Line_list(p9,1),2)==0.8)))
            bound_node2(p11) = [Line_list(p9,1)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,2),1)==1) && ((NL(Line_list(p9,2),2)==0.8)))
            bound_node2(p11) = [Line_list(p9,2)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,1),1)==1) && ((NL(Line_list(p9,1),2)==0.2)))
            bound_node2(p11) = [Line_list(p9,1)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,2),1)==1) && ((NL(Line_list(p9,2),2)==0.2)))
            bound_node2(p11) = [Line_list(p9,2)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,1),1)==1) && ((NL(Line_list(p9,1),2)==0.5))) %|| ((NL(Line_list(p9,2),1)==1) && ((NL(Line_list(p9,2),2)==0.5)))
            bound_node2(p11) = [Line_list(p9,1)];
            p11=p11+1;
            % elseif ((NL(Line_list(p9,2),1)==1) && ((NL(Line_list(p9,2),2)==0.5))) %|| ((NL(Line_list(p9,2),1)==1) && ((NL(Line_list(p9,2),2)==0.5)))
            %     bound_node2(p11) = [Line_list(p9,2)];
            %     p11=p11+1;
        elseif ((NL(Line_list(p9,1),1)==0.8) && ((NL(Line_list(p9,1),2)==0)))
            bound_node2(p11) = [Line_list(p9,1)];
            p11=p11+1;
        elseif ((NL(Line_list(p9,2),1)==0.8) && ((NL(Line_list(p9,2),2)==0)))
            bound_node2(p11) = [Line_list(p9,2)];
            p11=p11+1;
        end
        p27 = p27+1;
    end
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
% map=bsxfun(@eq, EL(Leftidx),Left_indices(:));
bcon = unique(Conduit_Surf);
no_el = size(Conduit_Surf,1);
no_nodes = length(bcon); bcon2 = sort(bcon);
Conduit_Surf_pseudo = zeros(no_nodes,2);

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
        dir_term_3rd(loc_idx,loc_idx,i) = -k*(Ly/6)*[2,1;1,2];   %3rd term
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
        diff_term(:,:,i) =-(1/(2*A(i))).*(b(i,:)'*b(i,:)+c(i,:)'*c(i,:));     %1st term
        adv_term(:,:,i) = -(1/3)*Pe(1)*[b(i,:);b(i,:);b(i,:)]+Pe(2)*[c(i,:);c(i,:);c(i,:)];     %2nd term
    end

end

%% Place the 2D FDM Model as a function in here
%% Assemble matrices and vectors []
for i = 1:size(EL,1)
    for i2 = 1:size(EL,2)
        for i3 = 1:size(EL,2)
            nno = EL(i,:);      % Return the node numbers in the element
            vN_temp = dir_term_3rd(i2,i3,i);
            diff_temp = diff_term(i2,i3,i);
            adv_temp = adv_term(i2,i3,i);
            K(nno(i2),nno(i3)) = K(nno(i2),nno(i3))+vN_temp+diff_temp+adv_temp;
            M(nno(i2),nno(i3)) = M(nno(i2),nno(i3))+m_el(i2,i3);
        end
        vN_temp_5th = von_neum_cond_5th(i2,1,i);
        dir_temp_4th = direchlet_cond_4th(i2,1,i);
        U(nno(i2),1)  = U(nno(i2),1) + dir_temp_4th + vN_temp_5th;
    end
end
%%
%
%   % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 2;
%   % set the seconds per image
% % open the video writer
%   open(writerObj);

t_real = 0;
for t = 2:no_t_steps
    %     if t ==1
    %          K_inv = pinv(K);
    %          C = K_inv*U;
    %     else
    t_real = t_real+dt;
    Peter = pinv(M-dt*(K));
    C_ot(:,t) = Peter*(M*C_ot(:,(t-1))+dt*(-U));
    %      M_inv = pinv(M);
    %      C_ot(:,t) = M_inv*((M*C_ot(:,t-1)+dt*(U-(K)*C_ot(:,t-1))));
    peter = figure(3)
    patch('Faces',EL,'Vertices',NL,'CData',C_ot(:,t),'FaceColor','interp');
    axis equal
    axis([0 1 0 1])
    hold on
    for i = 1:size(conduit_coordinates1,1)
        plot([conduit_coordinates1(i) conduit_coordinates3(i)],[conduit_coordinates2(i) conduit_coordinates4(i)],'LineWidth',1.8,'color','k')
        hold on
    end
    hold off
    P_title = sprintf('Concentration Plot t = %f',t_real)
    title(P_title);
    %      F(i) = getframe(gcf);
    drawnow
    %     end
    toc
end
end

