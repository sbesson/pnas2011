%% coleochaete_growth
% Treat the growth and division of Coleochaete
%%
%% Description
% Creates a cellNetwork object formed of a divided circle which
% grows and divides. In this program, we only consider the division of the
% central triangular cell.
%
%% Syntax   
% coleochaete_growth
%
%% Inputs
% none
%
%% Outputs
% 
%% Example
% >>  coleochaete_growth
%
%% Author
% Sebastien% Besson
% email address: sbesson@oeb.harvard.edu
% February 2008; Last revision: October 17, 2008

% Clean the existing variables and figures
clear
clf
close all

% Disable the debugging messages
warning off all

% Choose if a movie must be saved
save_movie = 0;
save_image=0;
count=0;
if(save_movie)
    aviobj = avifile('ColeochaeteGrowth.avi');
    frame_count = 1;
end

% Create new figure
f1 = figure('Position',[100 100 600 600],'Color',[1 1 1]);
axes('Position',[0 0 1 1]);
N=20;
exprate = 2^(1/(2*N));

V = [-1 0; 
    0 -1;
    1 0;
    0 1;];
E = [1 2 pi/4 1;
    2 3 pi/4 1;
    3 4 pi/4 1;
    4 1 pi/4 1;];
C = {[1 2 3 4];};
tissue = cellNetwork(V,E,C);
S0 = abs(cellArea(tissue,1));
h = plot(tissue,'-k');
axis([-21 21 -21 21]);
axis off
if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
if (save_image)
    count=count+1;
    I = getframe(gcf);
    imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
end

%while cellArea(tissue,1) <2*S0
for i=1:N
    tissue = expand(tissue,exprate,exprate,1);
    h = replot(tissue,h,'-k');
    drawnow;
    if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
    if (save_image)
        count=count+1;
        I = getframe(gcf);
        imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
    end
end

tissue.v=[tissue.v;0 0];
tissue.e = [tissue.e;1 5 0 2; 5 3 0 2];
tissue.c ={[1 2 -6 -5];[5 6 3 4]};
h = replot(tissue,h,'-k');
if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
if (save_image)
    count=count+1;
    I = getframe(gcf);
    imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
end

%while cellArea(tissue,1) <2*S0
for i=1:N
    tissue = expand(tissue,exprate,exprate,1);
    h = replot(tissue,h,'-k');
    drawnow;
    %waitforbuttonpress
    if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
    if (save_image)
        count=count+1;
        I = getframe(gcf);
        imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
    end
end

tissue.e = [tissue.e;2 5 0 3; 5 4 0 3];
tissue.c ={[1 7 -5];[2 -6 -7]; [3 -8 6];[8 4 5]};
h = replot(tissue,h,'-k');
if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end

%while cellArea(tissue,1) <2*S0
for i=1:N
    tissue = expand(tissue,exprate,exprate,1);
    h = replot(tissue,h,'-k');
    drawnow;
    %waitforbuttonpress
    if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
    if (save_image)
        count=count+1;
        I = getframe(gcf);
        imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
    end
end

planes = divisionplanes(tissue,1:4);
tissue =divide(tissue,[planes([1 4],1); planes([2 3],2)]);
h = replot(tissue,h,'-k');
if (save_movie), aviobj = addframe(aviobj,getframe(f1)); end
if (save_image)
    count=count+1;
    I = getframe(gcf);
    imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
end

while length(tissue.c)<348
    % Growth step
    tissue = expand(tissue,exprate,exprate,1);

    cells = find(abs(cellArea(tissue)) >2*S0);
    if ~isempty(cells)
        
        planes = divisionplanes(tissue,cells);
        tissue = divide(tissue,planes(:,1));     
        h = replot(tissue,h,'-k');
        %axis tight equal off;
        drawnow;
        %button = questdlg('Keep on growing?','Question','Yes','No','Yes');
        %si = all(button=='Yes') 
    end
    
    h = replot(tissue,h,'-k');
    %axis tight equal off;
    drawnow;
    % Add frame to movie
    if(save_movie), aviobj = addframe(aviobj,getframe(f1));end 
    if (save_image)
        count=count+1;
        I = getframe(gcf);
        imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
    end
end
if(save_movie), aviobj = close(aviobj); end 
if (save_image)
    count=count+1;
    I = getframe(gcf);
    imwrite(I.cdata,[ 'Movie/movie' num2str(count) '.png']);
end
plot2svg('Coleochaete-FinalPattern.svg',gcf);
return

% % 4 quadrant cells
V = [-pi/4 -pi/4; 
    0 0;
    pi/4 pi/4;
    -pi/4 pi/4;
    pi/4 -pi/4;];
E = [1 2 0 2;
    2 3 0 2;
    3 4 pi/4 1;
    4 1 pi/4 1;
    1 5 pi/4 1;
    5 3 pi/4 1;
    4 2 0 3;
    2 5 0 3;];
C = {[1 -7 4];
    [2 3 7];
    [5 -8 -1];
    [6 -2 8]};

% Original motif
% x=.3;
% V = [-cos(pi/4) -sin(pi/4); 
%     -x*cos(pi/4) -x*sin(pi/4); 
%     x*cos(pi/4) -x*sin(pi/4); 
% 	cos(pi/4) -sin(pi/4); 
%     cos(pi/4) sin(pi/4);
%     -cos(pi/4) sin(pi/4);];
% E = [1 2 0 1;
%     2 3 0 1;
%     3 4 0 1;
%     3 5 0 1;
%     2 6 0 1;
%     1 4 pi/4 1;
%     4 5 pi/4 1;
%     5 6 pi/4 1;
%     6 1 pi/4 1;];
% C = {[1 2 3 -6 ];
%     [3 7 -4];
%     [2 4 8 -5];
%     [1 5 9]};

% V = [-1 0; 
%      1 0;];
%  E = [1 2 0 1;
%      2 1 pi/2 1;
%      1 2 pi/2 1;];
%  C = {[1 2]; [3 -1]};
%  
tissue = cellNetwork(V,E,C);
S0 = abs(2*cellArea(tissue,1));

h = plot(tissue);
axis tight equal off;
% planes(1,:) = divisionplanes(tissue,1,.45);
% planes(2,:) = divisionplanes(tissue,2,.55);
planes = divisionplanes(tissue,1:4);
tissue =divide(tissue,[planes([1 4],1); planes([2 3],2)]);
%tissue = merge(tissue,[1 6]);
 %h = replot(tissue,h,'-k');
n_divisions = 2;
i = 1;

% Division cycles
while length(tissue.c)<348
    % Growth step
    tissue = expand(tissue,1.01,1.01,1);

    cells = find(abs(cellArea(tissue)) >2*S0);
    if ~isempty(cells)
        
        planes = divisionplanes(tissue,cells);
        tissue = divide(tissue,planes(:,1));     
        h = replot(tissue,h,'-k');
        %axis tight equal off;
        drawnow;
        %button = questdlg('Keep on growing?','Question','Yes','No','Yes');
        %si = all(button=='Yes') 
    end
    
    h = replot(tissue,h,'-k');
    %axis tight equal off;
    drawnow;
    % Add frame to movie
    if(save_movie), aviobj = addframe(aviobj,getframe(f1));end 
end

% Close the movie
if(save_movie), aviobj = close(aviobj); end 