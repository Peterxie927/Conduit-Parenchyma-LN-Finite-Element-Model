function [Ux Uy] = gradient_solver(Pressure_sol,NL,EL,NoN,Ly,width_x,width_y)
%% Parameters % Use cgs unit system for all the parameters.
mu = 1E-02; % [(dyne s)/cm^2]  Lymph viscosity
K_dimensionless =  [2.5E-08 2.5E-08]; % [cm^2] Permeability tensor
rou = 1000; % [kg/m^3] Fluid density
g = 9.81; % [m^2/s] Gravitational constant
% K_cond = (mu/K_dimensionless)*rou*g; % [] Pseudo hydraulic conductivity
f = 1000; % [cm^2/s] mass flow source term
dt = 1e-3; % Time-step [s]
total_time = 0.1; % Total time for PDE [s]
no_t_steps = total_time/dt; % Number of time-steps [ ]
% nx = 20; % Number of Nodes x-direction [ ]
% ny = 20; % Number of Nodes y-direction [ ]
% NoN = (nx+1)*(ny+1); % Total number of nodes 2D
L = 1e-2; % Characteristic Length
p9 = 1; p10 = 1; p11 = 1; p12 = 1;
q = 0; % Flux on Left Boundary (Von Neumann condition)


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
pres_term = zeros(3,3);
vel_term = zeros(3,1);

debug_1 = 0; debug_2 = 0; debug_3 = 0; debug_4 = 0; debug_5 = 0;
Kx = zeros(1,3);
Ky = zeros(1,3);
P = zeros(NoN,1);
pres_tempx = zeros(3,1); pres_tempy = zeros(3,1);
Ux = zeros(size(EL,1),1);
Uy = zeros(size(EL,1),1);

for i=1:size(EL,1)  % solve Weak Form PDE for each element
    pres_term_x(:,:) = (1/(2*A(i)))*(K_dimensionless(1)*[b(i,:);b(i,:);b(i,:)]);    %2nd term
    pres_term_y(:,:) = (1/(2*A(i)))*(K_dimensionless(2)*[c(i,:);c(i,:);c(i,:)]);
    for i2 = 1:3
            nno = EL(i,:);      % Return the node numbers in the element
            pres_tempx(i2) = pres_term_x(1,i2);
            pres_tempy(i2) = pres_term_y(1,i2);
            Kx(1,i2) = pres_tempx(i2);
            Ky(1,i2) = pres_tempy(i2);
            P_temp = [Pressure_sol(nno(1)) Pressure_sol(nno(2)) Pressure_sol(nno(3))];
    end
    Ux(i,1) = Kx*(P_temp)';
    Uy(i,1) = Ky*(P_temp)';
end
end



