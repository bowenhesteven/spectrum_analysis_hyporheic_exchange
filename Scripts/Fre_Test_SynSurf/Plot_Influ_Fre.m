%% Plot the Influence of the # fitted frequencies on approximating the real surface with different complexity
N_WL=[2,4,8,16,32,64];
plot(N_WL,mean_error_matrix1,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
hold on;
plot(N_WL,mean_error_matrix2,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','r',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
hold on;
plot(N_WL,mean_error_matrix4,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','m',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
hold on;
plot(N_WL,mean_error_matrix8,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','g',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
hold on;
plot(N_WL,mean_error_matrix16,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','y',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
hold on;
plot(N_WL,mean_error_matrix32,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',60)
hold on;
plot(N_WL,mean_error_matrix_real(1:6),'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','c',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',60)
hold on;


ax = gca;
ax.FontSize = 18; 
xlabel('Number of fitted frequencies','FontSize',20)
ylabel('Mean error','FontSize',20)
legend('N=1 Real Surface','N=2 Real Surface','N=4 Real Surface','N=8 Real Surface','N=16 Real Surface','N=32 Real Surface','Section of real surface')