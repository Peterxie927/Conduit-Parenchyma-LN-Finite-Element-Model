function [] = plot_mesh(NL,EL)
%UNTITLED Summary of this function goes here
figure (2)
trimesh(EL,NL(:,1),NL(:,2), 'Linewidth',2, 'Color',[0 0.4470 0.7410])
hold on;
plot (NL(:,1), NL(:,2), "r*");
nodes = [1:length(NL(:,1))]' ;  % nodes 
ele = [1:size(EL,1)]' ; % elements 
text(NL(:,1),NL(:,2),num2str(nodes),'color','b') 
% text(mean(NL(:,1)(tri),2),mean(y(tri),2),num2str(ele),'color','k')
set(gcf,'color','w')
xlabel('x [m]')
ylabel('y [m]')
hold on
end

