clear
close all
clc

% Load Dionaea quadrant results
load('./dionaea_results.mat');


f1 = figure('Position',[20 500 1000 500]);
msize =  10;

drawTriangularPhaseDiagram;
hold on

theta1 = arrayfun(@(x) x.modelshape.a1, dionaea_results);
theta1(theta1<0) = -theta1(theta1<0);
theta2 = arrayfun(@(x) x.modelshape.a2, dionaea_results);
theta2(theta2<0) = -theta2(theta2<0);
plot(theta1([dionaea_results.divisionmode] == 1),...
    theta2([dionaea_results.divisionmode] == 1),...
    'or','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize)
plot(theta1([dionaea_results.divisionmode] == 2),...
    theta2([dionaea_results.divisionmode] == 2),...
    'or','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize)
plot(theta1([dionaea_results.divisionmode] == 3),...
    theta2([dionaea_results.divisionmode] == 3),...
    'or','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize)
