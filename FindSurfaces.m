function [Left_Surf,Left_connectivity,Left_indices,Right_Surf,Right_connectivity,Right_indices,Bottom_Surf, Bottom_connectivity, Bottom_indices,Top_Surf, Top_connectivity,Top_indices] = FindSurfaces(NL,EL,width_x,width_y)
% FIND_SURFACES Find the surfaces present in the mesh of the finite element
% model
p1 = 1; p2 = 1; p3 = 1; p4 = 1;
for i = 1:size(NL,1)
  if NL(i,1) == 0
    Left_Surf(p1,:) = NL(i,:);
    Left_connectivity(p1,:) = EL(i,:);
    Left_indices(p1) = i;
    p1 = p1+1;
  elseif NL(i,1) == width_x
    Right_Surf(p2,:) = NL(i,:);
    Right_connectivity(p2,:) = EL(i,:);
    Right_indices(p2) = i;
    p2 = p2+1;
  elseif NL(i,2) == 0
    Bottom_Surf(p3,:) = NL(i,:);
    Bottom_connectivity(p3,:) = EL(i,:);
    Bottom_indices(p3) = i;
    p3 = p3+1;
  elseif NL(i,2) == width_y
    Top_Surf(p4,:) = NL(i,:);
    Top_connectivity(p4,:) = EL(i,:);
    Top_indices(p4) = i;
    p4 = p4+1;
  disp(i)
end
end

