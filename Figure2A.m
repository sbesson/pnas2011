%Figure2A.m
clear
close all
clc

% Load Coleochaete results
load('./coleochaeteresults.mat');


f1 = figure('Position',[20 500 1000 500]);
msize =  10;

drawColeochaetePhaseDiagram;
hold on

theta = arrayfun(@(x) x.modelshape.theta, coleochaeteresults);
h = arrayfun(@(x) x.modelshape.h, coleochaeteresults);
plot(theta([coleochaeteresults.divisionmode] == 1),...
    h([coleochaeteresults.divisionmode] == 1),...
    'or','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',msize)
plot(theta([coleochaeteresults.divisionmode] == 2),...
    h([coleochaeteresults.divisionmode] == 2),...
    'or','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',msize)
plot(theta([coleochaeteresults.divisionmode] == 3),...
    h([coleochaeteresults.divisionmode] == 3),...
    'or','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',msize)
