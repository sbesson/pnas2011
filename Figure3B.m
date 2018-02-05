clear
close all
clc

% Load Dionaea quadrant results
load('./quadrantcells.mat');


f1 = figure('Position',[20 500 1000 500]);
msize =  10;

drawQuadrantPhaseDiagram;
hold on

theta1 = arrayfun(@(x) x.modelshape.a1, quadrantcells);
theta2 = arrayfun(@(x) x.modelshape.a2, quadrantcells);
plot(theta1([quadrantcells.divisionmode] == 1),...
    theta2([quadrantcells.divisionmode] == 1),...
    'or','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize)
plot(theta1([quadrantcells.divisionmode] == 2),...
    theta2([quadrantcells.divisionmode] == 2),...
    'or','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize)
plot(theta1([quadrantcells.divisionmode] == 3),...
    theta2([quadrantcells.divisionmode] == 3),...
    'or','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize)
