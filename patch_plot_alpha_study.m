clear all
close all
clc

load('Conduit_med_profiles.mat')
t_real = 0;
figure (1)
for t = 9999
        t_real = t*dt;
        P_title = sprintf('t = %f',t_real+dt);
        title(P_title);
        patch('Faces',EL,'Vertices',NL,'CData',C_ot(1:NoN,t),'FaceColor','interp','Edgecolor','none');
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
        colormap turbo

        
end