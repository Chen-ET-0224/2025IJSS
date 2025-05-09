function Sphere_indenter_VdW 
%w0                           %dimensionless indentation depth
K=100;                        
r_indenter=200;               %dimensionless indenter radius
tau_array=[0,25,250];
col=parula(length(tau_array)+1);
nu=0.3;                       %poissons ratio
NUM=1500;                     
infinity=500;                 
conarea=5;                    %dimensionless contact area

for i=1:length(tau_array)
tau=tau_array(i);                                  %dimensionless pretension
for j=1:1                     
%------------first guess of w0
guess1=0.1;
w0=guess1;
solinit = bvpinit(logspace(log10(conarea),log10(infinity),NUM),@ex1init);
options = bvpset('stats','on');
sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
for s=1:3
    sol = bvp4c(@ex5ode,@ex5bc,sol,options);
end
r = sol.x;
f = sol.y;
deri=f(3,:)-f(2,:)./r./r;
deri1=deri(1)-2/r_indenter;
%------------second guess of w0
guess2=0.11;
w0=guess2;
solinit=sol;
options = bvpset('stats','on');
sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
for s=1:3
    sol = bvp4c(@ex5ode,@ex5bc,sol,options);
end
r = sol.x;
f = sol.y;
deri=f(3,:)-f(2,:)./r./r;
deri2=deri(1)-2/r_indenter;

while(abs(deri2)>=0.00001)
    temp=guess2-0.5*(deri2/(deri2-deri1))*(guess2-guess1);
    guess1=guess2;
    guess2=temp;
%calculation after iteration----------------------------------------------
deri1=deri2;
w0=guess2;
solinit=sol;
options = bvpset('stats','on');
sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
for s=1:3
    sol = bvp4c(@ex5ode,@ex5bc,sol,options);
end
r = sol.x;
f = sol.y;
inter=linspace(0.001,conarea,300);
sphere=-w0+(inter.*inter)./(r_indenter);
force=abs(trapz(r,-2*pi*vdWW(f(1,:)).*r))+abs(trapz(inter,-2*pi*vdWW(sphere).*inter));
deri=f(3,:)-f(2,:)./r./r;
deri2=deri(1)-2/r_indenter;
end

x=linspace(0.001,conarea,100);
y=-w0+x.^2./r_indenter;
plot(x,y/w0,'Color',col(i,:),'LineWidth',1,'HandleVisibility','off');
% plot(x,y/w0,'Color',col(i,:),'LineWidth',1.5,'HandleVisibility','off');
hold on 
plot(-x,y/w0,'Color',col(i,:),'LineWidth',1,'HandleVisibility','off');

plot(r,f(1,:)/w0,'Color',col(i,:),'LineWidth',1,'HandleVisibility','off');
hold on;
plot(-r,f(1,:)/w0,'Color',col(i,:),'LineWidth',1);
axis([-15 15 -1.4 0.8])
conarea=5
w0


%---------------------------------Equation solver
rs=r_indenter;
lp=(tau+sqrt(tau^2-4))/2;
ln=(tau-sqrt(tau^2-4))/2;
lp=lp^0.5;
ln=ln^0.5;

%------------first guess
guess1=0.1;
eta=guess1;
k0p=besselk(0,lp*eta);
k0n=besselk(0,ln*eta);
k1p=besselk(1,lp*eta);
k1n=besselk(1,ln*eta);
k2p=besselk(2,lp*eta);
k2n=besselk(2,ln*eta);
c2=-(-2*k0n*eta-k1n*ln*eta*eta+w0*k1n*ln*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
c4=-(2*k0p*eta+k1p*lp*eta*eta-w0*k1p*lp*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
V=0.5*c2*lp*lp*(k0p+k2p)+0.5*c4*ln*ln*(k0n+k2n);
deri1=V-2/rs;
%------------second guess
guess2=0.105;
eta=guess2;
k0p=besselk(0,lp*eta);
k0n=besselk(0,ln*eta);
k1p=besselk(1,lp*eta);
k1n=besselk(1,ln*eta);
k2p=besselk(2,lp*eta);
k2n=besselk(2,ln*eta);
c2=-(-2*k0n*eta-k1n*ln*eta*eta+w0*k1n*ln*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
c4=-(2*k0p*eta+k1p*lp*eta*eta-w0*k1p*lp*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
V=0.5*c2*lp*lp*(k0p+k2p)+0.5*c4*ln*ln*(k0n+k2n);
deri2=V-2/rs;



while(abs(deri2)>=0.001)
    temp=guess2-0.1*(deri2/(deri2-deri1))*(guess2-guess1);
    guess1=guess2;
    guess2=temp;
    %calculation after iteration----------------------------------------------
    deri1=deri2;
    eta=guess2;
    k0p=besselk(0,lp*eta);
    k0n=besselk(0,ln*eta);
    k1p=besselk(1,lp*eta);
    k1n=besselk(1,ln*eta);
    k2p=besselk(2,lp*eta);
    k2n=besselk(2,ln*eta);
    c2=-(-2*k0n*eta-k1n*ln*eta*eta+w0*k1n*ln*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
    c4=-(2*k0p*eta+k1p*lp*eta*eta-w0*k1p*lp*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
    V=0.5*c2*lp*lp*(k0p+k2p)+0.5*c4*ln*ln*(k0n+k2n);
    deri2=V-2/rs;
end
eta 
w0
x1=linspace(0.01,eta,100);
y1=-w0+x1.^2./rs;
plot(x1,y1/w0,'k--','HandleVisibility','off','LineWidth',1);
hold on
plot(-x1,y1/w0,'k--','HandleVisibility','off','LineWidth',1);
x2=linspace(eta,infinity,2000);
lp=(tau+sqrt(tau^2-4))/2;
ln=(tau-sqrt(tau^2-4))/2;
lp=lp^0.5;
ln=ln^0.5;
k0p=besselk(0,lp*eta);
k0n=besselk(0,ln*eta);
k1p=besselk(1,lp*eta);
k1n=besselk(1,ln*eta);
k2p=besselk(2,lp*eta);
k2n=besselk(2,ln*eta);
c2=-(-2*k0n*eta-k1n*ln*eta*eta+w0*k1n*ln*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
c4=-(2*k0p*eta+k1p*lp*eta*eta-w0*k1p*lp*rs)/(rs*(k0p*k1n*ln-k0n*k1p*lp));
y2=c2.*besselk(0,lp.*x2)+c4.*besselk(0,ln.*x2);
plot(x2,y2/w0,'k--','HandleVisibility','off','LineWidth',1);
plot(-x2,y2/w0,'k--','HandleVisibility','off','LineWidth',1);
conarea=5;         
ax = gca;
    ax.TickLabelInterpreter = 'latex';
    grid on
    ax.MinorGridAlpha = 0.1;
    ax.GridAlpha = 0.1;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickLength = [0.02 0.02];
    set(ax,'fontsize',10)
    set(ax,'LineWidth',0.5)
    xlabel('$R=r/l$','interpreter','latex','fontsize',20);
    ylabel('$W=w/\delta$','interpreter','latex','fontsize',20);
end
end
legend({'$\tau=0$','$\tau=25$','$\tau=250$'},'Interpreter','latex');
% --------------------------------------------------------------------------
function vdW=vdWW(r)
    % vdW=(8/3).*(1./(r+1).^3-1./(r+1).^9);
    vdW=(1/6).*(1./(r+1).^3-1./(r+1).^9);
end
function dfdr = ex5ode(r,f)   
%f(1)=w
%f(2)=r*w'
%f(3)=(r*w')'/r
%f(4)=((r*w')'/r)'*r
%f(5)=psi
%f(6)=psi'
%f(7)=psi''  

dfdr = [ f(2)./r
          r.*f(3)
          f(4)./r
          (-r.*((1/6).*(1./(f(1)+1).^3-1./(f(1)+1).^9))+f(6).*(f(3)-f(2)./r./r)+f(2).*f(7)./r)
          f(6)
          f(7)
          -f(7)./r+f(6)./r./r-K.*f(2).*f(2)./r./r./r 
        ];
end
% -------------------------------------------------------------------------
function res = ex5bc(f0,finf)
res = [f0(1)+w0-(conarea*conarea)/(r_indenter)                                %pre-defined indentation depth
       f0(2)/conarea-(2*conarea)/r_indenter                                   %Slope continuity condition
       f0(5)                                                                  %zero stress function


       % f0(7)-nu.*f0(6)./conarea                                             %Non-slipping conditions

       f0(7)-f0(6)./conarea+K*(conarea*conarea)/(r_indenter*r_indenter)       %perfect-slip conditions

       finf(1)                                                                
       finf(2)./infinity                                                      %far-field conditions
       % finf(7)-nu*finf(6)/infinity                                         
       finf(7)-nu.*finf(6)./infinity-(1-nu)*tau  
      ];
end
function y = ex1init(r)
lp=1i;
lm=-1i;
ypred=-2*w0/log(lm/lp)*(besselk(0,lp^0.5*r)-besselk(0,lm^0.5*r));
    
y = [  ypred
           0
           0
           0
           0
           0
           0
     ];
    
end
end