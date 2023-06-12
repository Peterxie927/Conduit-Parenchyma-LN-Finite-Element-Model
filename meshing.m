function [N_List,E_List] = meshing(width_x,width_y,nx,ny,element_type);

Dim = 2; % Number of dimensions in the problem
boundary = [0 0; width_x 0; 0 width_y; width_x width_y];  % Store corners of the domain

NoN = (nx+1)*(ny+1); % Total number of nodes
NoE = nx*ny; % Total number of Elements
NoPoints = 4;

%% Define Nodes
N_List = zeros(NoN,Dim);
a = width_x/nx; % Increment in x-direction
b = width_y/ny; % Increment in y-direction

n=1;

for i = 1:ny+1
    for j =1:nx+1
        N_List(n,1) = boundary(1,1) + (j-1)*a;  % j is columns in node list. x-coordinates
        N_List(n,2) = boundary(1,2) + (i-1)*b;  % i is rows in node list. y-xoordinates
        n = n+1;
    end
end

%% Define Elements
E_List = zeros(NoE,NoPoints);

for i = 1:ny
    for j = 1:nx
        if j==1
            E_List((i-1)*nx+j,1) = (i-1)*(nx+1)+j;
            E_List((i-1)*nx+j,2) = E_List((i-1)*nx+j,1)+1;
            E_List((i-1)*nx+j,4) = E_List((i-1)*nx+j,1) + (nx+1);
            E_List((i-1)*nx+j,3) = E_List((i-1)*nx+j,4)+1;
        else
            E_List((i-1)*nx+j,1) = E_List((i-1)*nx+j-1,2);
            E_List((i-1)*nx+j,4) = E_List((i-1)*nx+j-1,3);
            E_List((i-1)*nx+j,2) = E_List((i-1)*nx+j,1) + 1;
            E_List((i-1)*nx+j,3) = E_List((i-1)*nx+j,4)+1;
        end
    end
end

if isequal(element_type,'triangle')
   N_Points_tri = 3;
   EList_new = zeros(2*NoE,N_Points_tri);

   for i = 1:NoE
       EList_new(2*(i-1)+1,1) = E_List(i,1);
       EList_new(2*(i-1)+1,2) = E_List(i,2);
       EList_new(2*(i-1)+1,3) = E_List(i,3);

       EList_new(2*(i-1)+2,1) = E_List(i,1);
       EList_new(2*(i-1)+2,2) = E_List(i,3);
       EList_new(2*(i-1)+2,3) = E_List(i,4);
   end

   E_List = EList_new;
end

