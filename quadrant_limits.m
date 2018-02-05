function quadrant_limits
% Clean the existing variables and figures
clear
clf
close all
clc
% Disable the debugging messages
warning off all

dtheta =pi/64;
angle = 90;
savefile = 'quadrant_limits.mat';
if exist(savefile,'file') ==2, load(savefile); end

eqline1 = [pi/4 pi/4; pi/3 pi/3];
eqline2 = [pi/3 pi/3; 5*pi/16 3*pi/8];
eqline3 = [pi/3 pi/3; 5*pi/8 3*pi/16];
save(savefile,'eqline1','eqline2','eqline3');

levels = [.08 .16];
l_curves = struct();
for i= 1:length(levels)
    l = levels(i);
    l_curves(i).level = l;
    display(['Level ' num2str(l)]);
    smin = fminsearch(@(s) (Ldiff(s,[-1 1 0])-l)^2 + (Ldiff(s,[0 1 -1]))^2,[pi/3 pi/3],optimset('MaxIter',50));
    
    theta1c = smin(1);
    theta2c = smin(2);
   
    display('Curve1');
    x1 = pi/4:dtheta:theta1c; 
    y1 = arrayfun(@(theta1) fzero(@(theta2) Ldiff([theta1 theta2],[-1 1 0])-l,theta1),x1);

    y2 = 3*pi/8:-dtheta:theta2c; 
    x2 = arrayfun(@(theta2) fzero(@(theta1) Ldiff([theta1 theta2],[-1 0 1])-l,pi/2-theta2/2),y2);

    l_curves(i).curve1 =  [x1 theta1c x2(end:-1:1);y1 theta2c y2(end:-1:1)]';
%     plot(x1,y1,'or');
%     plot(x2,y2,'or');
    
    smin = fminsearch(@(s)  (Ldiff(s,[1 -1 0])-l)^2 + (Ldiff(s,[0 -1 1])-l)^2,[pi/3 pi/3],optimset('MaxIter',50));
    theta1c = smin(1);
    theta2c = smin(2);
    display('Curve2');
    x1 = pi/4:dtheta:theta1c;
    %y1 = pi/8:dtheta:theta2c; 
    y1 = arrayfun(@(theta1) fzero(@(theta2) Ldiff([theta1 theta2],[1 -1 0])-l,theta1),x1);

    x2 = 5*pi/8:-dtheta:theta1c; 
    y2 = arrayfun(@(theta1) fzero(@(theta2) Ldiff([theta1 theta2],[0 -1 1])-l,pi/2-theta1/2),x2);

    l_curves(i).curve2 =  [x1 theta1c x2(end:-1:1);y1 theta2c y2(end:-1:1)]';
    
    smin = fminsearch(@(s) (Ldiff(s,[1 0 -1])-l)^2 + (Ldiff(s,[0 1 -1])-l)^2,[pi/3 pi/3],optimset('MaxIter',50));
    theta1c = smin(1);
    theta2c = smin(2);
   display('Curve3');
    y1 = 3*pi/8:-dtheta:theta2c;
    x1 = arrayfun(@(theta2) fzero(@(theta1) Ldiff([theta1 theta2],[1 0 -1])-l,pi/2-theta2/2+pi/8),y1);

    x2 = 5*pi/8:-dtheta:theta1c; 
    y2 = arrayfun(@(theta1) fzero(@(theta2) Ldiff([theta1 theta2],[0 1 -1])-l,pi/2-theta1/2),x2);

    l_curves(i).curve3 =  [x1 theta1c x2(end:-1:1);y1 theta2c y2(end:-1:1)]';
    
    save(savefile,'eqline1','eqline2','eqline3','l_curves');
end


    function L = divisionlengths(s)
        T = tri2cell(triangle(s(1),s(2)),angle*pi/180);
        area = abs(cellArea(T));
        M = divisionplanes(T);
        L = NaN(1,3);
        i1 = find(([M.i] == 2) & ([M.j] == 3));
        if ~isempty(i1)
            L(1) = M(i1).walllength/sqrt(area);
        end
        i2 = find(([M.i] == 1) & ([M.j] == 3));
        if ~isempty(i2)
            L(2) = M(i2).walllength/sqrt(area);
        end
        i3 = find(([M.i] == 1) & ([M.j] == 2));
        if ~isempty(i3)
            L(3) = M(i3).walllength/sqrt(area);
        end
    end

%     function delta =findtriplepoint(s,L,l)
%         X = L.*divisionlengths(s);
%         delta = (X(3)-X(1)-l)^2+(X(2)-X(1)-l)^2;
%     end
        
    function delta = Ldiff(s,L)
        X = L.*divisionlengths(s);
        delta = sum(X(L~=0));
    end

end