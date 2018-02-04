function  [S,L,L2,T] = drawColeochaetePhaseDiagram(facealpha)

if nargin==0, facealpha =.33; end


taucr = fzero(@(tau) findtaucr(tau),pi/4);
function delta=findtaucr(tau)
    h2 = (pi/2-tau)*(tan(tau));
    theta2 = 2*((pi/2-tau)*(tan(tau))^2+tau-tan(tau))/(h2.*(2-h2));
    delta=theta2-sqrt(h2^2./(1-h2+(h2^2)/2));
end
hcr = (pi/2-taucr)*(tan(taucr));

% Draw first line of equilibrium
h1=0:0.0001:hcr;
theta1 = sqrt(h1.^2./(1-h1+(h1.^2)/2));

% Draw second line of equilibrium
tau = taucr:.002:pi/2;
h2 = (pi/2-tau).*(tan(tau));
theta2 = 2*((pi/2-tau).*(tan(tau)).^2+tau-tan(tau))./(h2.*(2-h2));

% Draw third line of equilibrium
taumax = fzero( @(tau) sqrt(2)*((pi/2-tau).*(tan(tau)).^2+tau-tan(tau))-(pi/2-tau)*tan(tau),pi/4);
tau = [taucr:.002:taumax taumax];
At =@(tau) (pi/2-tau).*(tan(tau)).^2+tau-tan(tau);
theta3 = At(tau)/2 + sqrt(At(tau).^2+4*(pi/2-tau).^2.*tan(tau).^2)/2;
h3 = 1- sqrt(2*((pi/2-tau).*tan(tau)./theta3).^2-1);
h3(end) =1;

hold on

S(1) = fill([theta1 theta2 pi pi], [h1 h2 1 0],'b','FaceAlpha',facealpha);
S(2) = fill([theta1 theta3 0 0], [h1 h3 1 0],'r','FaceAlpha',facealpha);
S(3) = fill([theta3(end:-1:1) theta2 pi], [h3(end:-1:1) h2 1],'g','FaceAlpha',facealpha);

L(1) = plot(theta1,h1,'-k','LineWidth',3);
L(2) = plot(theta2,h2,'-k','LineWidth',3);
L(3) = plot(theta3,h3,'-k','LineWidth',3);

levels = [.08 .16];

A = @(theta,h) theta/2*h*(2-h);
Lp = @(theta,h) theta*sqrt(1/2+1/2*(1-h)^2);
La = @(theta,h) h;
Ld = @(tau) (pi/2-tau).*tan(tau);
At = @(tau) (pi/2-tau).*(tan(tau)).^2+tau-tan(tau);
TAU = [pi/4 pi/4+pi/42];
for i= 1:length(levels)
    l=levels(i);

    % Lines 1
    theta_12 = @(tau) At(tau)/2+At(tau)/2.*(1+4*(Ld(tau)./At(tau)-l).^2).^.5;
    h_12 = @(tau) 1- (1-2*At(tau)./theta_12(tau)).^.5;
    tau0 = fzero(@(tau) La(theta_12(tau),h_12(tau))-Lp(theta_12(tau),h_12(tau))-l*sqrt(At(tau)),pi/8);
    
    h11 = 0:0.002:h_12(tau0);
    theta11 = arrayfun(@(h) fzero(@(theta) La(theta,h)-Lp(theta,h)-l*sqrt(A(theta,h)),h),h11);
   
    theta12 = @(tau) At(tau)/2+At(tau)/2.*(1+4*(Ld(tau)./At(tau)-l).^2).^.5;
    h12 = @(tau) 1- (1-2*At(tau)./theta12(tau)).^.5;
    taumax = fzero(@(tau) 2*At(tau)-theta12(tau),tau0);
    tau = tau0:.002:taumax;
    L2(i,1) = plot(theta11,h11,'--k','LineWidth',1);
    L2(i,2) = plot([theta12(tau) theta12(taumax)],[h12(tau) 1],'--k','LineWidth',1);

    % Lines 2
    theta_21 = @(tau) At(tau)/2+At(tau)/2.*(1+4*(Ld(tau)./At(tau)+l).^2).^.5;
    h_21 = @(tau) 1- (1-2*At(tau)./theta_21(tau)).^.5;
    tau0 = fzero(@(tau) La(theta_21(tau),h_21(tau))-Ld(tau)-l*sqrt(At(tau)),TAU(i));
   
    theta21 = @(tau) At(tau)/2+At(tau)/2.*(1+4*(Ld(tau)./At(tau)+l).^2).^.5;
    h21 = @(tau) 1- (1-2*At(tau)./theta21(tau)).^.5;
    h22 = @(tau) Ld(tau)+l*sqrt(At(tau));
    theta22 = @(tau) 2*At(tau)./(h22(tau).*(2-h22(tau)));
    
    taumax = fzero(@(tau) 2*At(tau)-theta21(tau),tau0);
    tau =  tau0:.002:taumax;
    L2(i,3) = plot([theta21(tau) theta21(taumax)],[h21(tau) 1],'--k','LineWidth',1);

    taumax = fzero(@(tau) 1-Ld(tau)-l*sqrt(At(tau)),tau0);
    tau =  [tau0:.002:taumax taumax];
    L2(i,4) = plot(theta22(tau),h22(tau),'--k','LineWidth',1);
    
    % Lines 3
    theta_31 = @(tau) At(tau)/2+At(tau)/2.*(1+4*(Ld(tau)./At(tau)).^2).^.5;
    h_31 = @(tau) 1- (1-2*At(tau)./theta_31(tau)).^.5;
    tau0 = fzero(@(tau) Ld(tau)- La(theta_31(tau),h_31(tau))-l*sqrt(At(tau)),pi/8);
    
    h31 = 0:0.002:h_31(tau0);
    theta31 = arrayfun(@(h) fzero(@(theta) Lp(theta,h)-La(theta,h)-l*sqrt(A(theta,h)),h),h31);    
    h32 = @(tau) Ld(tau)-l*sqrt(At(tau));
    theta32 = @(tau) 2*At(tau)./(h32(tau).*(2-h32(tau)));
    
    taumax = fzero(@(tau) theta32(tau)-pi,pi/4);
    tau = tau0:.002:taumax;
	L2(i,5) = plot(theta31,h31,'--k','LineWidth',1);
    L2(i,6) = plot(theta32(tau),h32(tau),'--k','LineWidth',1);
    
    theta0 = 3*pi/4+pi/16;
    tau0 = fzero(@(tau) theta32(tau)-theta0,pi/4);
    %h_text = interp1(theta32(tau0),h32(tau0),theta_text);
    T(i) = text(theta32(tau0)-.1,h32(tau0)+.04,['\delta_{12}=' num2str(l)],'FontSize',18,'FontName','Arial');
end

axis([0 pi 0 1]);
axis square
set(gca,'FontName','Helvetica','FontSize',20)
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'},'XTick',0:pi/4:pi)
set(gca,'YTick',0:.2:1);
ylabel('H/R','FontSize',24,'FontName','Helvetica')
xlabel('\theta','FontSize',24,'FontName','Helvetica')
end