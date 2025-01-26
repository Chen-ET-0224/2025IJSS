function cylindrical_test_vdw3
infinity = 500;   %A dimensionless length used to represent infinity
conarea=1e-6;     %Coordinate origin (avoid singularity);
r_indenter=10;
p=1e9;      %Stiffness of imaginary spring
fac=1;
eta=r_indenter;
K1=100;
K2=100;
K3=100;
nu=0.3;           %POISSON'S RATIO
tau_array=logspace(log10(0.001),log10(100),100);
ktau=[];
tautau=[];
kout=[];
col=parula(length(tau_array)+1);
c=1;
vv=1;   %Coefficient selection variable (1: point model 0: spring model)



for i=1:length(tau_array)
    ddelta=0.9;        %Dimensionless distance between the indenter and the horizontal reference plane of the membrane
    delta=[];
    fn=[];
    tau=tau_array(i)
    for n=1:1
        t=1;
        force1=0.1;
        deltaF1_original=0.1; deltaF2_original=0.1; deltaF3_original=0.1;
        guess=[1-ddelta (1-ddelta)/2 (1-ddelta)/4];

        while (abs(deltaF1_original)>force1*0.001||abs(deltaF2_original)>force1*0.001||abs(deltaF3_original)>force1*0.001)
            lengthguess=length(guess);
            if t==1
                w0=guess(1);  w1=guess(2); w2=guess(3);

                % mesh1=logspace(log10(conarea),log10(r_indenter/5),300);
                % mesh2=logspace(log10(r_indenter/5),log10(infinity),100);
                % subset2 = mesh2(2:length(mesh2));
                % mesh=[mesh1,subset2];
                mesh1=linspace(conarea,r_indenter,300);
                mesh2=logspace(log10(r_indenter),log10(infinity),150);
                subset2 = mesh2(2:length(mesh2));
                mesh=[mesh1,subset2];
                solinit = bvpinit(mesh,@ex1init);
                options = bvpset('stats','on');

                if n==1&&i==1
                    sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
                end

                for s=1:2
                    sol = bvp4c(@ex5ode,@ex5bc,sol,options);
                end
                r = sol.x;   f = sol.y;
                force1=trapz(r,-2*pi*(           vdW1(r,(ddelta-f(1,:)-1))                ).*r);
                force2=trapz(r,-2*pi*vdW2(f(1,:)-f(8,:)).*r);
                force3=trapz(r,-2*pi*vdW2(f(8,:)-f(15,:)).*r);
                force4=trapz(r,-2*pi*vdW2(f(15,:)).*r);
                deltaF1_original=force1-force2;
                deltaF2_original=force2-force3;
                deltaF3_original=force3-force4;
                deltaF=[deltaF1_original deltaF2_original deltaF3_original];
            else
                deltaF=[deltaF1_original deltaF2_original deltaF3_original];
            end




            for v=1:lengthguess        %write Jacobian

                guessn=guess;
                guessn(v)=guessn(v)-0.0000001;
                w0=guessn(1);  w1=guessn(2);  w2=guessn(3);
                solinit = sol;
                options = bvpset('stats','on');
                sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
                for s=1:2
                    sol = bvp4c(@ex5ode,@ex5bc,sol,options);
                    options = bvpset('stats','on');
                end
                r = sol.x;    f = sol.y;
                force1=trapz(r,-2*pi*(           vdW1(r,(ddelta-f(1,:)-1))                ).*r);
                force2=trapz(r,-2*pi*vdW2(f(1,:)-f(8,:)).*r);
                force3=trapz(r,-2*pi*vdW2(f(8,:)-f(15,:)).*r);
                force4=trapz(r,-2*pi*vdW2(f(15,:)).*r);
                deltaF1=force1-force2;
                deltaF2=force2-force3;
                deltaF3=force3-force4;
                jacobian(1,v)=(deltaF1-deltaF1_original)/(guessn(v)-guess(v));
                jacobian(2,v)=(deltaF2-deltaF2_original)/(guessn(v)-guess(v));
                jacobian(3,v)=(deltaF3-deltaF3_original)/(guessn(v)-guess(v));
            end
            guess_reverse=guess'-fac*jacobian\deltaF';
            guesskk=guess_reverse';
            guess=guesskk;

            w0=guess(1);  w1=guess(2);  w2=guess(3);
            solinit=sol;
            options = bvpset('stats','on');
            sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
            for s=1:2
                sol = bvp4c(@ex5ode,@ex5bc,sol,options);
            end
            r = sol.x;
            f = sol.y;
            % figure(1)
            % plot(r,f(1,:));
            force1=trapz(r,-2*pi*(           vdW1(r,(ddelta-f(1,:)-1))                ).*r);
            force2=trapz(r,-2*pi*vdW2(f(1,:)-f(8,:)).*r);
            force3=trapz(r,-2*pi*vdW2(f(8,:)-f(15,:)).*r);
            force4=trapz(r,-2*pi*vdW2(f(15,:)).*r);
            deltaF1_original=force1-force2;
            deltaF2_original=force2-force3;
            deltaF3_original=force3-force4;
            t=t+1;
        end
        delta(n)=1-ddelta;
        fn(n)=force4;
        % if i==1
        %     figure(1)
        % plot(r,f(1,:));
        % hold on;
        % plot(r,f(8,:)-1);
        % plot(r,f(15,:)-2);
        % end
        %
        % if i==2
        %     figure(2)
        % plot(r,f(1,:));
        % hold on;
        % plot(r,f(8,:)-1);
        % plot(r,f(15,:)-2);
        % end
        %
        % if i==3
        %     figure(3)
        % plot(r,f(1,:));
        % hold on;
        % plot(r,f(8,:)-1);
        % plot(r,f(15,:)-2);
        % end
        %
        %     ax = gca;
        %     ax.XLim = [0 50];
        % title('Indentation');
        % xlabel('r');
        % ddelta=ddelta-0.03;
    end

    ktau(i)=(fn(1))/(delta(1));
    tautau(i)=tau;

    if vv==1
        algha=0.4110;
    else
        algha=1/3;
    end
    % algha=1/3;
    etas=algha^(1/4)*eta;
    taus=tau*algha^(-1/2);
    lp=(taus+sqrt(taus^2-4*c))/2;
    ln=(taus-sqrt(taus^2-4*c))/2;
    lp=lp^0.5;
    ln=ln^0.5;
    k0p=besselk(0,lp*etas);
    k0n=besselk(0,ln*etas);
    k1p=besselk(1,lp*etas);
    k1n=besselk(1,ln*etas);
    c2d=-k1n*ln/(k0p*k1n*ln-k0n*k1p*lp);
    c4d=-k1p*lp/(-k0p*k1n*ln+k0n*k1p*lp);
    kout(i)=sqrt(algha)*abs(2*pi*c*etas*(c2d/lp*k1p+c4d/ln*k1n));
    if isnan(kout(i))==1||isinf(kout(i))==1
        c2dd=-k1n.*ln./(k1n.*ln-k0n.*lp);
        c4dd=-lp./(-k1n.*ln+k0n.*lp);
        kout(i)=sqrt(algha)*abs(2*pi*c*etas*(c2dd./lp+c4dd./ln.*k1n));
    end
    %     algha=0.4012;
    % etas=algha^(1/4)*eta;
    % taus=tau*algha^(-1/2);
    %     k_limit(i)=sqrt(algha)*2*pi*sqrt(taus.^2 - 4*c)./atanh(sqrt(1 - 4*c./taus./taus));







    % figure(4)
    % loglog(delta,fn,'Color',col(i,:),'LineWidth',1.5);
    % hold on;
end


figure(3)
loglog(tautau,ktau,'b')
hold on
loglog(tautau,kout+sqrt(algha)*pi*etas^2,'r')
% loglog(tautau,k_limit,'g')
kks=kout+sqrt(algha)*pi*etas^2;
ks=ktau;
errorL = max(abs((kks - ks) ./ ks))


% Functions--------------------------------------------------------------------------
    function vdW=vdW1(r,w)
        vdW = zeros(size(r));


        for o = 1:length(r)
            if r(o) > r_indenter
                vdW(o) = 0;
            else
                vdW(o) = p * w(o);
            end
        end
    end


    function vdWW=vdW2(r)
        % vdWW=(1/6).*(1./(r+1).^3-1./(r+1).^9);
        vdWW=r;
    end


    function dfdr = ex5ode(r,f)

        %f(1)=w
        %f(2)=r*w'
        %f(3)=(r*w')'/r
        %f(4)=((r*w')'/r)'*r
        %f(5)=psi
        %f(6)=psi'
        %f(7)=psi''

        %f(8)=w2
        %f(9)=r*w2'
        %f(10)=(r*w2')'/r
        %f(11)=((r*w2')'/r)'*r
        %f(12)=psi2
        %f(13)=psi2'
        %f(14)=psi2''

        %f(15)=w3
        %f(16)=r*w3'
        %f(17)=(r*w3')'/r
        %f(18)=((r*w3')'/r)'*r
        %f(19)=psi3
        %f(20)=psi3'
        %f(21)=psi3''
        dfdr = [  f(2)./r
            r.*f(3)
            f(4)./r
            -r.*(        vdW2(f(1)-f(8))-vdW1(r,(ddelta-f(1)-1)) )+f(6).*(f(3)-f(2)./r./r)+f(2).*f(7)./r
            f(6)
            f(7)
            -f(7)./r+f(6)./r./r-K1.*f(2).*f(2)./r./r./r

            f(9)./r
            r.*f(10)
            f(11)./r
            -r.*(        vdW2(f(8)-f(15))-vdW2(f(1)-f(8)) )+f(13).*(f(10)-f(9)./r./r)+f(9).*f(14)./r
            f(13)
            f(14)
            -f(14)./r+f(13)./r./r-K2.*f(9).*f(9)./r./r./r


            f(16)./r
            r.*f(17)
            f(18)./r
            -r.*(        vdW2(f(15))-vdW2(f(8)-f(15)) )+f(20).*(f(17)-f(16)./r./r)+f(16).*f(21)./r
            f(20)
            f(21)
            -f(21)./r+f(20)./r./r-K3.*f(16).*f(16)./r./r./r
            ];
    end
% -------------------------------------------------------------------------
    function res = ex5bc(f0,finf)
        res = [f0(1)+w0
            f0(2)/conarea
            f0(5)
            f0(7)-nu.*f0(6)./conarea
            finf(1)
            finf(2)/infinity
            finf(7)-nu.*finf(6)./infinity-(1-nu)*tau

            f0(8)+w1
            f0(9)/conarea
            f0(12)
            f0(14)-nu.*f0(13)./conarea
            finf(8)
            finf(9)/infinity
            finf(14)-nu.*finf(13)./infinity-(1-nu)*tau


            f0(15)+w2
            f0(16)/conarea
            f0(19)
            f0(21)-nu.*f0(20)./conarea
            finf(15)
            finf(16)/infinity
            finf(21)-nu.*finf(20)./infinity-(1-nu)*tau
            ];
    end
    function y = ex1init(r)
        ypred1=-w0;
        ypred2=-w1;
        ypred3=-w2;

        y = zeros(size(r));


        for e = 1:length(r)
            if r(e) > r_indenter
                lp=sqrt(c)*1i;
                ln=-sqrt(c)*1i;
                lp=sqrt(lp);
                ln=sqrt(ln);
                ypred1(e)=-2*w0/log(ln/lp)*(besselk(0,lp^0.5*r)-besselk(0,ln^0.5*r));
                ypred2(e)=-2*w1/log(ln/lp)*(besselk(0,lp^0.5*r)-besselk(0,ln^0.5*r));
                ypred3(e)=-2*w2/log(ln/lp)*(besselk(0,lp^0.5*r)-besselk(0,ln^0.5*r));
            else
                ypred1(e)=-w0;
                ypred2(e)=-w1;
                ypred3(e)=-w2;
            end
        end


        y = [        ypred1
            0
            0
            0
            0
            0
            0
            ypred2
            0
            0
            0
            0
            0
            0
            ypred3
            0
            0
            0
            0
            0
            0
            ];
    end
end