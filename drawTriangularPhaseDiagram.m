function  [S,L,L2] = drawTriangularPhaseDiagram(facealpha)

if nargin==0, facealpha =.33; end

load('triangular_limits.mat');

hold on

S(1) = fill([eqline1(:,1)' eqline3(:,1)' pi/2 eqline1(1,1)],[eqline1(:,2)' eqline3(:,2)' pi/8 eqline1(1,2)],'b','FaceAlpha',facealpha);
S(2) = fill([eqline1(:,1)' eqline2(:,1)' pi/8 pi/8 eqline1(1,1)],[eqline1(:,2)' eqline2(:,2)' pi/2 pi/8 eqline1(1,2)],'r','FaceAlpha',facealpha);
S(3) = fill([eqline3(:,1)' pi/2 eqline2(end:-1:1,1)'],[eqline3(:,2)' pi/2 eqline2(end:-1:1,2)'],'g','FaceAlpha',facealpha);

L(1) = plot(eqline1(:,1),eqline1(:,2),'-k','LineWidth',3);
L(2) = plot(eqline2(:,1),eqline2(:,2),'-k','LineWidth',3);
L(3) = plot(eqline3(:,1),eqline3(:,2),'-k','LineWidth',3);
    
for i=1:length(l_curves) 
    L2(i,1) = plot(l_curves(i).curve1(:,1),l_curves(i).curve1(:,2),'--k','LineWidth',1);
    L2(i,2) = plot(l_curves(i).curve2(:,1),l_curves(i).curve2(:,2),'--k','LineWidth',1);
    L2(i,3) = plot(l_curves(i).curve3(:,1),l_curves(i).curve3(:,2),'--k','LineWidth',1);
    T(i) = text(pi/8+.1,i*pi/8+.04,['\delta_{12}=' num2str(l_curves(i).level)],...
        'Rotation',45,'FontSize',18,'FontName','Arial');
end


% Add graphical components
axis([pi/8 pi/2 pi/8 pi/2]);
axis square
set(gca,'FontName','Helvetica','FontSize',20)
set(gca,'XTickLabel',{'\pi/8','\pi/4','3\pi/8','\pi/2'},'XTick',pi/8:pi/8:pi/2)
set(gca,'YTickLabel',{'\pi/8','\pi/4','3\pi/8','\pi/2'},'YTick',pi/8:pi/8:pi/2)
ylabel('\theta2','FontSize',24)
xlabel('\theta1','FontSize',24)

end