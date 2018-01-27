%% quadrantcell_movie

% Clean the existing variables and figures
clear
clc
close all

% Disable the debugging messages
warning off all

% Initialize figure
fig=figure('Position',[100 100 800 450],'Color',[1 1 1]);
axes('OuterPosition',[0 0 1/2 1]);

T = tri2cell(triangle(pi/2,pi/4+pi/48),pi/2);
load('quadrantcell-solutions.mat');
hold on

% I = min([planes(:).i],[planes(:).j]);
% J = max([planes(:).i],[planes(:).j]);
s = [planes.arcposition];
l = [planes.length];
I = [planes(:).i];
J =[planes(:).j];
s1 = edgelength(T,1)/sum(edgelength(T));
s2 = sum(edgelength(T,1:2))/sum(edgelength(T));

%index = find([planes(:).i]==[planes(:).j]);
%plot([planes(index).arcposition],[planes(index).length],'ok');
hold on;

index11a = (I==1 & J==I & s<.2);
h11a =plot(s(index11a),l(index11a),'-k');
index11b = (I==1 & J==I & s>.2);
h11b = plot(s(index11b),l(index11b),'-k');
index22a = (I==2 & J==2 & s<.6);
h22a = plot(s(index22a),l(index22a),'-k');
index22b = (I==2 & J==2 & s>.6);
h22b = plot(s(index22b),l(index22b),'-k');

index12 = (I==1 & J==2);
h12 = plot(s(index12),l(index12),'-','Color',coloredge(3));
index21 = (I==2 & J==1);
h21 =plot(s(index21),l(index21),'-','Color',coloredge(3));

index13 = (I==1 & J==3);
h13 =plot(s(index13),l(index13),'-','Color',coloredge(2));
index31 = (I==3 & J==1);
h31 = plot(s(index31),l(index31),'-','Color',coloredge(2));

index23 = (I==2 & J==3);
h23 =plot(s(index23),l(index23),'-','Color',coloredge(1));
index32 = (I==3 & J==2);
h32 = plot(s(index32),l(index32),'-','Color',coloredge(1));

% Add graphical components
plot([s1 s1],[.5 2],'--k');
plot([s2 s2],[.5 2],'--k');
axis([0 1 .5 2]);

set(gca,'FontName','Arial','FontSize',16)
ylabel('Length','FontSize',20,'FontName','Arial')
xlabel('Arc position','FontSize',20,'FontName','Arial')

set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%title('Energy landscape of a quadrant cell','FontSize',22)

% Add event functions to the axes and draws a point
axlim=axis;
currentpoint = plot([0 0],[axlim(3) axlim(4)],'--','Linewidth',3,'Color',[0 0 0],'Markersize',8);

axes('OuterPosition',[1/2 0 1/2 1]);
cellplot = plot(T,'-k','Linewidth',4);
for i =1:6
    divisionplot(i) = plot(cellNetwork(),'-k','Linewidth',2);
end
axis tight equal off;

set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

uniques=unique([planes(:).arcposition]);
count = 1;
for si=uniques
   
    %cellplot = replot(T,cellplot,'-k');

    index = find(s==si);
    for i=1:length(index)
        % Draw division plane 1
        if planes(index(i)).i == planes(index(i)).j
            color = [0 0 0];
        else
            color = coloredge(setdiff([1 2 3],...
                [planes(index(i)).i planes(index(i)).j]));
        end
        V = [planes(index(i)).Q; planes(index(i)).P];
        E = [1 2 planes(index(i)).angle 1];
        divisionplot(i) = replot(cellNetwork(V,E,{}),divisionplot(i),'Color',color);
    end

    for i=length(index)+1:length(divisionplot)
        divisionplot(i) = replot(cellNetwork(),divisionplot(i));
    end

    % Determine plane of division of shortest length
    [junk, imin] = min([planes(index).length]);
    set([divisionplot(:).edges],'Linewidth',2);
    set(divisionplot(imin).edges,'Linewidth',4);
    set(currentpoint,'XData',[si si],'HandleVisibility','on',...
            'EraseMode','normal','ButtonDownFcn',@movePoint_Start);
    axis equal off
    hold off
    set(h11a,'XData',s(index11a & s<=si),'YData',l(index11a & s<=si));
    set(h11b,'XData',s(index11b & s<=si),'YData',l(index11b & s<=si));
    set(h22a,'XData',s(index22a & s<=si),'YData',l(index22a & s<=si));
    set(h22b,'XData',s(index22b & s<=si),'YData',l(index22b & s<=si));
    set(h12,'XData',s(index12 & s<=si),'YData',l(index12 & s<=si));
    set(h21,'XData',s(index21 & s<=si),'YData',l(index21 & s<=si));
    set(h13,'XData',s(index13 & s<=si),'YData',l(index13 & s<=si));
    set(h31,'XData',s(index31 & s<=si),'YData',l(index31 & s<=si));
    set(h23,'XData',s(index23 & s<=si),'YData',l(index23 & s<=si));
    set(h32,'XData',s(index32 & s<=si),'YData',l(index32 & s<=si));
    drawnow;
    %waitforbuttonpress
    I = getframe(gcf);
%     imwrite(I.cdata, ['Images/Movie' num2str(count) '.jpg'], 'Quality', 80);
    count=count+1;

end
% Close the movie
