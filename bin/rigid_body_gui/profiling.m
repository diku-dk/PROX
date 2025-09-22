%close all;
%clear all;
%clc;







%red   = [0.7,0.1,0.1];
%green = [0.1,0.7,0.1];
%blue  = [0.1,0.1,0.7];
%f_type = 'Times';
%f_size = 20;


%frames = [1:length(Ekin)];
%Emech = Epot + Ekin;


%figure(100);
%e1 = plot(frames, Ekin,'r-'); 
%hold on; 
%e2 = plot(frames, Epot, 'g-'); 
%e3 = plot(frames, Emech, 'b-'); 
%hold off; 
%title('Kinetic and potential energy of total system', 'FontSize', f_size, 'FontName', f_type); 
%legend([e1,e2,e3], {'Kinetic', 'Potential','Mechanical'}, 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Energy (Joules)', 'FontSize', f_size, 'FontName', f_type)
%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type)
%print(gcf,'-depsc2','energy');
%print(gcf,'-dpng','energy');


%figure(200);
%plot(abs(average_penetration)*100);
%title('Penetrations', 'FontSize', f_size, 'FontName', f_type); 
%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Penetration depth (cm)', 'FontSize', f_size, 'FontName', f_type)
%print(gcf,'-depsc2','penetration');
%print(gcf,'-dpng','penetration');


%figure(300);
%plot(contacts);
%title('Contacts', 'FontSize', f_size, 'FontName', f_type); 
%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Number of contacts', 'FontSize', f_size, 'FontName', f_type)
%print(gcf,'-depsc2','contacts');
%print(gcf,'-dpng','contacts');


%figure(400);
%plot(contacts(contacts>0) ,solver, '.'); 
%title('Time used by solver vs number of contacts', 'FontSize', f_size, 'FontName', f_type); 
%xlabel('Number of contacts', 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Solver (ms)', 'FontSize', f_size, 'FontName', f_type)
%print(gcf,'-depsc2','solver_time_vs_contacts');
%print(gcf,'-dpng','solver_time_vs_contacts');


%figure(500);
%plot(frames(contacts>0), solver); 
%title('Time used by solver per frame', 'FontSize', f_size, 'FontName', f_type)
%xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Solver (ms)', 'FontSize', f_size, 'FontName', f_type)
%print(gcf,'-depsc2','solver_time_per_frame');
%print(gcf,'-dpng','solver_time_per_frame');


%figure(600);
%semilogy(cell2mat(convergence(1)));
%hold on
%for frame = 2:numel(convergence); 
%  semilogy(cell2mat(convergence(frame)));
%end
%title('Convergence', 'FontSize', f_size, 'FontName',f_type );
%xlabel('Solver iteration', 'FontSize', f_size, 'FontName', f_type); 
%ylabel('Natural merit function', 'FontSize', f_size, 'FontName', f_type);
%print(gcf,'-depsc2','convergence');
%print(gcf,'-dpng','convergence');


%figure(700);
%plot(cell2mat(rfactor(1)));
%hold on
%for frame = 2:numel(rfactor);
%  plot(cell2mat(rfactor(frame)));
%end
%title('R-factor development', 'FontSize', f_size, 'FontName',f_type );
%xlabel('Solver iteration', 'FontSize', f_size, 'FontName', f_type); %ylabel('R-factor', 'FontSize', f_size, 'FontName', f_type);%print(gcf,'-depsc2','rfactors');
%print(gcf,'-dpng','rfactors');


