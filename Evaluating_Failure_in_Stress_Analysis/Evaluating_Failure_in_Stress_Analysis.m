clc; close all; clear all;
%% Seo and Seligson 
%% Evaluating_Failure_in_Stress_Analysis_in_Infinite_Plate_with_Circular_Hole 

%% initialize 
double maxth;
double maxr;
double maxt; 
 
%radius of hole
a=1;
 
length=6*a; 
width=6*a;
 
%number of pieces each unit is divided into
div=100;
 
%resizing length in regards to division of pieces
l=length*div;
 
%degrees of surface which will be evaluated
d=360;
 
%initialize arrays
sigth=zeros([d l-a*div]);
sigr= zeros([d l-a*div]);
tau= zeros([d l-a*div]);
 
%% polar matrix for stress
 
%evaluate stress across surface (polar)
for theta=0:d
    for x=a*div:l
        r=x/div; 
   %need to index
        sigth(theta+1, x-(a*div-1))=((1+(a/r)^2)-(1+3*a^4/r^4)*cos(deg2rad(2*theta)))*0.5;
        sigr(theta+1, x-(a*div-1))=((1-(a/r)^2)+(1-4*(a/r)^2+3*a^4/r^4)*cos(deg2rad(2*theta)))/2;
        tau(theta+1, x-(a*div-1))=-(1+2*(a/r)^2-3*(a/r)^4)*sin(deg2rad(2*theta))/2;
    end
end
 
sige=sqrt(abs(sigr).^2 - abs(sigr.*sigth) + abs(sigth).^2 + 3.*abs(tau));
%% max stress
maxth=max(max(abs(sigth)));
 
maxr=max(max(abs(sigr)));
 
maxt=max(max(abs(tau)));
 
maxe=max(max(abs(sige)));
 
%% covert to cartesian
 
[r,t]=meshgrid(a:1/div:length, 0:pi/180:deg2rad(d));
x=r.*cos(t);
y=r.*sin(t);
 
pt1=linspace(a,length,(l-(a*div-1)));
 
%% plot
 
%question 1
figure
plot(pt1, sige(1,:), pt1,sige(31,:),pt1,sige(61,:),pt1,sige(91,:)) 
legend('stress at 0\circ','stress at 30\circ','stress at 60\circ','stress at 90\circ')
%only r/a when a=1
xlabel('r/a');
ylabel('\sigma_e/s');
xlim([a 4*a])
 
%question 2
figure 
contourf(x,y,sige);
xlim([-(4*a) (4*a)])
ylim([-(4*a) (4*a)]);
xlabel('x/a');
ylabel('y/a');
 
%question 3
disp('maximum (Von Mises stress)/s:'); disp(maxe);



