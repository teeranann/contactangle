function [Data]=EllipticFitVolnotilt(edgeL,edgeR,baseL,baseR,nupline)
% Function that fits the dropboundaries edgeL and edgeR to two seperate
% ellipses. Then the contact angle is measured as the angle of the ellipse
% at the crossing between the fitted ellipse and the baseline defined by the
% two points baseL and baseR.
% 
% 
% 
% Cut point below n up line
edgeLN.x=edgeL.x(edgeL.y < baseL(2)-nupline);
edgeLN.y=edgeL.y(edgeL.y < baseL(2)-nupline);
edgeRN.x=edgeR.x(edgeR.y < baseL(2)-nupline);
edgeRN.y=edgeR.y(edgeR.y < baseL(2)-nupline);
%---------------------------------------------------------------
%        Rename & restructure input variables
%---------------------------------------------------------------

left=[edgeLN.x,edgeLN.y];
right=[edgeRN.x,edgeRN.y];
% left=[edgeL.x,edgeL.y]; 
% right=[edgeR.x,edgeR.y];


x0L=baseL(1);
y0L=baseL(2);

x0R=baseR(1);
y0R=baseR(2);

%---------------------------------------------------------------
%         Fit ellipses to the data above the baseline
%---------------------------------------------------------------

indexL=find_index(left,baseR,baseL);
indexR=find_index(right,baseR,baseL);

%assign empty structure to ellipse in case fits fails.
efL.a=[];
efR.a=[];
% efL=fit_ellipsenotilt(left(1:indexL,1),left(1:indexL,2));
% efR=fit_ellipsenotilt(right(1:indexR,1),right(1:indexR,2));
efL=fit_ellipsenotilt(left(:,1),left(:,2));
efR=fit_ellipsenotilt(right(:,1),right(:,2));

if ~isempty(efL.a) && ~isempty(efR.a)
% Find tripple line as crossing between baseline and ellipse
[O1,~]=EllipseLineIntersection(efL,[x0L,y0L],[x0R,y0R]);
[~,O2]=EllipseLineIntersection(efR,[x0L,y0L],[x0R,y0R]);

% Calculate left angle between drop and x axis and put data into output structure

alpha_indexL=atan2(O1(2)-efL.Y0_in,O1(1)-efL.X0_in);
t_indexL=alpha2t(alpha_indexL,efL);
CAL=calculatecontactangle(efL,t_indexL,1);

Data.CA_flat{1}=CAL;
Data.TLL=[O1(1),O1(2)];
Data.t0{1}=t_indexL;
Data.ellipse{1}=efL;

% Calculate right angle between drop and x axis and put data into output structure
alpha_indexR=atan2(O2(2)-efR.Y0_in,O2(1)-efR.X0_in);
t_indexR=alpha2t(alpha_indexR,efR);
CAR=calculatecontactangle(efR,t_indexR,2);

Data.CA_flat{2}=CAR;
Data.TLR=[O2(1),O2(2)];
Data.t0{2}=t_indexR;
Data.ellipse{2}=efR;

%Calculate true contact angle by compensating for the tilt of the baseline
tilt=-atan((O2(2)-O1(2))/(O2(1)-O1(1)))*180/pi; %minus sign is due to coordinates being upside down
Data.tilt=tilt;
    
Data.CAL=Data.CA_flat{1}-tilt;
Data.CAR=Data.CA_flat{2}+tilt;

Data.trace{1}=left;
Data.trace{2}=right;

else % if no ellipse was fitted fill in empty data structure
    for lr=1:2
Data.CA_flat{lr}=NaN;
Data.basepoint{lr}=[NaN,NaN];
Data.t0{lr}=NaN;
Data.basepoint{lr}(1)=NaN;
Data.basepoint{lr}(2)=NaN;
Data.CA_flat{lr}=NaN;
Data.ellipse{lr}.X0_in=NaN;
Data.ellipse{lr}.Y0_in=NaN;
Data.ellipse{lr}.a=NaN;
Data.ellipse{lr}.b=NaN;
Data.ellipse{lr}.phi=NaN;
    Data.tilt=NaN;
    Data.CA{1}=NaN;
    Data.CA{2}=NaN;
    end
end





function [O1,O2]=EllipseLineIntersection(ellipse,P1,P2)
a=ellipse.a;
b=ellipse.b;
x0=ellipse.X0_in;
y0=ellipse.Y0_in;
phi=-ellipse.phi;

myfun=@(x) root2d(x,x0,y0,phi,a,b,P1,P2) ;

t1=alpha2t(atan2(P1(2)-y0,P1(1)-x0),ellipse);
xstart = [t1,0];
options = optimoptions(@fsolve,'Display','off');
x1 = fsolve(myfun,xstart,options);

t2=alpha2t(atan2(P2(2)-y0,P2(1)-x0),ellipse);
xstart = [t2,1];
x2 = fsolve(myfun,xstart,options);

[Xe1,Ye1]=ellipserim(ellipse,x1(1));
[Xe2,Ye2]=ellipserim(ellipse,x2(1));
O1=[Xe1,Ye1];
O2=[Xe2,Ye2];

function F = root2d(x,X0,Y0,phi,a,b,P1,P2)
t=x(1);
s=x(2);

F(1) = X0+a*cos(t)*cos(phi)-b*sin(t)*sin(phi)-P1(1)-s*(P2(1)-P1(1));
F(2) = Y0+a*cos(t)*sin(phi)+b*sin(t)*cos(phi)-P1(2)-s*(P2(2)-P1(2));

function CA=calculatecontactangle(ellipse,t0,lr)
a=ellipse.a;
b=ellipse.b;

phi=-ellipse.phi;

Xtm=(-a*sin(t0)*cos(phi)-b*cos(t0)*sin(phi));
Ytm=(-a*sin(t0)*sin(phi)+b*cos(t0)*cos(phi));

if lr==1
CA=-atan2(Ytm,Xtm)*180/pi;
elseif lr==2
CA=atan2(Ytm,Xtm)*180/pi;
end

function [Xe,Ye]=ellipserim(ellipse,t)

a=ellipse.a;
b=ellipse.b;
x0=ellipse.X0_in;
y0=ellipse.Y0_in;
phi=-ellipse.phi;

Xe=x0+cos(t)*cos(phi)*a-sin(t)*sin(phi)*b;
Ye=y0+cos(t)*sin(phi)*a+sin(t)*cos(phi)*b;

function t0=alpha2t(alpha,ellipse)
a=ellipse.a;
b=ellipse.b;
phi=-ellipse.phi;
t0=atan2(a*(sin(alpha)*cos(phi)-cos(alpha)*sin(phi)),b*(cos(alpha)*cos(phi)+sin(alpha)*sin(phi)));