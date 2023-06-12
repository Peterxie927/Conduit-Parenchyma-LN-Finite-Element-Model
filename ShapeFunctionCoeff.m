function [A,a,b,c] = Shape_function_coeff(NL,EL)
%SHAPE_FUNCTION_COEFF Summary of this function goes here

for i=1:size(EL,1)
    A(i) = (1/2).*abs(NL(EL(i,1),1)*((NL(EL(i,2),2))-NL(EL(i,3),2))+NL(EL(i,2),1)*((NL(EL(i,3),2))-NL(EL(i,1),2))+NL(EL(i,3),1)*((NL(EL(i,1),2))-NL(EL(i,2),2)));
end
% Find coefficient 'a' in shape function
for i = 1:size(EL,1)
    a(i,1)=NL(EL(i,2),1)*NL(EL(i,3),2)-NL(EL(i,3),1)*NL(EL(i,2),2);
    a(i,2)=NL(EL(i,3),1)*NL(EL(i,1),2)-NL(EL(i,1),1)*NL(EL(i,3),2);
    a(i,3)=NL(EL(i,1),1)*NL(EL(i,2),2)-NL(EL(i,2),1)*NL(EL(i,1),2);
end
% Find coefficient 'b' in shape function
for i = 1:size(EL,1)
    b(i,1)=NL(EL(i,2),2) - NL(EL(i,3),2);
    b(i,2)=NL(EL(i,3),2) - NL(EL(i,1),2);
    b(i,3)=NL(EL(i,1),2) - NL(EL(i,2),2);
end
% Find coefficient 'c' in shape function
for i = 1:size(EL,1)
    c(i,1)=NL(EL(i,3),1) - NL(EL(i,2),1);
    c(i,2)=NL(EL(i,1),1) - NL(EL(i,3),1);
    c(i,3)=NL(EL(i,2),1) - NL(EL(i,1),1);
end

% for i=1:1
%     N1 = 1/(2*A(i))*(a(i,1)+b(i,1)

end

