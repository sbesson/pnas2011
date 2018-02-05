%% coleochaete_phasediagram
% Plot the phase diagram for the division of Coleochaete cells
%%
%% Description
% Based on the theoretical approach of the division of four-sided
% Coleochaete cells, this program plots the phase diagram of division of
% these cells as a function of two geometrical parameters : the angle theta
% and the ratio of lengths H/R.
%
% The programm allows to select point on the phase diagram and the right
% panel draws the corresponding cell with the allowed planes of division.
% The plane of division with the smallest length is drawn with a thicker
% line.
%
%
%% Syntax   
% coleochaete_phasediagram
%
%% Inputs
% none
%
%% Outputs
% 
%% Example
% >>  coleochaete_phasediagram
%
%% Author
% Sebastien% Besson
% email address: sbesson@oeb.harvard.edu
% August 2008; Last revision: November 21, 2008

function coleochaete_phasediagram(f)

% Clean the existing variables and figures
clear
clc
close all
% Disable the debugging messages
warning off all

if nargin<1, f=1; end
    
theta0 = pi/4;
h0 = .2;

% Initialize figure
fig=figure('Position',[50 50 2000 1000],'Color',[1 1 1]);

axes('OuterPosition',[0 0 1/2 1]);
%[S,L,E,G] = drawPhasediagram(f);
%set(G(:),'Visible','off');
%set(S(:),'ButtonDownFcn',@drawCell);
drawColeochaetePhaseDiagram(f);
currentpoint = plot(theta0,h0,'o','MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[1 1 1],'Markersize',8,...
    'ButtonDownFcn',@movePoint_Start,'EraseMode','normal');

axes('OuterPosition',[1/2 0 1/2 1]);
[C Divisions] = coleochaeteCell(theta0,h0);
cellhandle = plot(C,'-k','Linewidth',3);
divisionhandle = plot(Divisions(1),'--','Linewidth',3);            
[junk index] = min(edgelength(Divisions));
set(divisionhandle.edges(index),'LineStyle','-','Linewidth',4);
axis equal off

    function movePoint_Start(src,eventdata)
        set(currentpoint(:),'EraseMode','xor');
        set(fig,'WindowButtonMotion',@drawCell,'WindowButtonUp',@movePoint_Finish)
    end

    function movePoint_Finish(src,eventdata)
        drawCell(src,eventdata)
        set(currentpoint(:),'EraseMode','normal');
        set(fig,'WindowButtonMotion','','WindowButtonUp','');
    end

    function [C D] = coleochaeteCell(theta,h)   
       
        % Create and draw cell network for the mother cell
        C = col2cell(coleochaete(h,theta));
        R = sqrt(2/(theta*h*(2-h))); % radius
        d = 4/3*R*sin(theta/2)/theta*(1+(1-h)^2/(2-h)); % shift to origin        
        rn=sqrt(.08/(2*pi)*theta*h*(2-h))*R; % nucleus radius
       
        % Soap bubble model
        V=[];
        E=[];
        if h<1
            % Draw division plane 1
            V = [R-d 0; R*(1-h)-d 0];
            E = [1 2 0 3];
        end

        % Draw division plane 2
        Rd = R*fzero(@(l) theta/2*l^2-1/2*theta*(1-h)^2-theta/4*h*(2-h),.5);
        if (Rd<=R)
            V = [V;
                Rd*cos(theta/2)-d -Rd*sin(theta/2);
                Rd*cos(theta/2)-d Rd*sin(theta/2);];
            E = [E;
                size(V,1)-1 size(V,1) theta/2 1];
        end
        
        % Determine position of division plane 3
        tau = fzero(@(tau) (pi/2-abs(tau))*tan(abs(tau))^2+...
            abs(tau)-tan(abs(tau))-theta/2*h*(2-h),pi/6);
        % Test existence of division plane 3
        if (tau<=theta) && (tan(tau)*(1-sin(tau))+1-cos(tau)<=h) 
            % Draw division plane 3
            Rd = R*cos(tau)-R*tan(tau)*(1-sin(tau));
            V = [V;
                Rd*cos(theta/2)-d -Rd*sin(theta/2);
                R*cos(-theta/2+tau)-d R*sin(-theta/2+tau)];
            E = [E; size(V,1)-1 size(V,1) -(pi/2-tau)/2 2];
        end
        D = cellNetwork(V,E,{});
    end

    function drawCell(src,event)        
        % Get position of the pointer
        P = get(gca,'CurrentPoint');
        A = get(gca,{'XLim','YLim'});
        % Normalize coordinates within the axis limits
        theta = min(max(P(1,1),A{1}(1)), A{1}(2));
        h = min(max(P(2,2),A{2}(1)), A{2}(2));
        % Update position of the current point
        set(currentpoint(:),'XData',theta,'YData',h);

        [C Divisions] = coleochaeteCell(theta,h);
        
        cellhandle = replot(C,cellhandle,'-k','Linewidth',3);
        divisionhandle = replot(Divisions,divisionhandle); 
        set([divisionhandle.edges],'LineStyle','--','Linewidth',2);
        
        [junk index] = sort(edgelength(Divisions));
        set(divisionhandle.edges(index(1)),'LineStyle','-','Linewidth',6);
        set(divisionhandle.edges(index(2)),'LineStyle','--','Linewidth',4);

    end
end