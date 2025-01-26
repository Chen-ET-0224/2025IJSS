function cylindrical_test_vdw
infinity = 500;   %A dimensionless length used to represent infinity
conarea=1e-6;     %Coordinate origin (avoid singularity);
r_indenter=0.001;
p=1e11;        %Stiffness of imaginary spring
eta=r_indenter;
K=100;
c=1;
nu=0.3;           %POISSON'S RATIO

tau_array=logspace(log10(0.001),log10(100),100);
ktau=[];
tautau=[];
kout_theo=[];

for t=1:length(tau_array)
    ddelta=0.95;        %Dimensionless distance between the indenter and the horizontal reference plane of the membrane
    delta=[];
    fn=[];
    tau=tau_array(t)




    for n=1:1
        if n==1
            guess1=0.1;
        else
            guess1=guess2;
        end
        w0=guess1;       %Estimated membrane center deflection
        % mesh1=logspace(log10(conarea),log10(r_indenter),300);
        mesh1=linspace(conarea,r_indenter,300);
        mesh2=logspace(log10(r_indenter),log10(infinity),100);
        subset2 = mesh2(2:length(mesh2));
        mesh=[mesh1,subset2];
        solinit = bvpinit(mesh,@ex1init);
        options = bvpset('stats','on');
        if n==1&&t==1
            sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
        end
        for s=1:3
            sol = bvp4c(@ex5ode,@ex5bc,sol,options);
        end
        r = sol.x;
        f = sol.y;



        % for l=length(r)
        %     if r(l)>=r_indenter
        %         index=l;
        %         break;
        %     end
        % end
        % rr=r(index:length(r));
        % ff=f(1,index:length(r))




        force1=trapz(r,-2*pi*(           vdW1(r,1e5*(-f(1,:)-1+ddelta))/1e5                ).*r);
        force2=trapz(r,-2*pi*vdW2(f(1,:)).*r);
        deltaF1=force1-force2;








        guess2=guess1-0.0002;
        w0=guess2;       %Estimated membrane center deflection
        solinit=sol;
        options = bvpset('stats','on');
        sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
        for s=1:3
            sol = bvp4c(@ex5ode,@ex5bc,sol,options);
        end
        r = sol.x;
        f = sol.y;
        force1=trapz(r,-2*pi*(           vdW1(r,1e5*(-f(1,:)-1+ddelta))/1e5                ).*r);
        force2=trapz(r,-2*pi*vdW2(f(1,:)).*r);
        deltaF2=force1-force2;
        while(abs(deltaF2)>=force1*0.00001)
            temp=guess2-1*(deltaF2/(deltaF2-deltaF1))*(guess2-guess1);
            guess1=guess2;
            guess2=temp;
            deltaF1=deltaF2;
            w0=guess2;
            solinit = sol;
            options = bvpset('stats','on');
            sol = bvp4c(@ex5ode,@ex5bc,solinit,options);
            for loop=1:3
                sol = bvp4c(@ex5ode,@ex5bc,sol,options);
            end
            r = sol.x;
            f = sol.y;
            force1=trapz(r,-2*pi*(           vdW1(r,1e5*(-f(1,:)-1+ddelta))/1e5                ).*r);
            force2=trapz(r,-2*pi*vdW2(f(1,:)).*r);
            deltaF2=force1-force2;
        end
        force=trapz(r,-2*pi*vdW2(f(1,:)).*r);

        delta(n)=1-ddelta;
        fn(n)=force1;


        figure(1)
        plot(r,f(1,:));
        hold on;
        axis([0 25 -1 0.5])
        title('Indentation');
        xlabel('r');
        ylabel('displacement');
        % ddelta=ddelta-0.02;
    end

    ktau(t)=(fn(1))/(delta(1));
    tautau(t)=tau;


    lp=(tau+sqrt(tau^2-4*c))/2;
    ln=(tau-sqrt(tau^2-4*c))/2;
    lp=lp^0.5;
    ln=ln^0.5;
    k0p=besselk(0,lp*eta);
    k0n=besselk(0,ln*eta);
    k1p=besselk(1,lp*eta);
    k1n=besselk(1,ln*eta);
    c2d=-k1n*ln/(k0p*k1n*ln-k0n*k1p*lp);
    c4d=-k1p*lp/(-k0p*k1n*ln+k0n*k1p*lp);
    kout(t)=abs(2*pi*c*eta*(c2d/lp*k1p+c4d/ln*k1n));

end
force1


figure(3)
loglog(tautau,ktau,'b','LineWidth',1.5)
hold on
loglog(tautau,kout+pi*r_indenter^2,'r')









% loglog(tautau,pi*r_indenter^2                +8.*tautau./tautau,'k--')






% FUNCTIONS--------------------------------------------------------------------------
% function vdW=vdW1(r,w)
%     if r>r_indenter
%         vdW=0;
%     else
%         vdW=p*w;
%     end
% end

    function vdW=vdW1(r,w)
        vdW = zeros(size(r));


        for v = 1:length(r)
            if r(v) > r_indenter
                vdW(v) = 0;
            else
                vdW(v) = p*w(v);
                % vdW(v) = (-0.5*p*(    r(v)/r_indenter    )+2*p)*w(v);
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
        dfdr = [  f(2)./r
            r.*f(3)
            f(4)./r
            -r.*(        vdW2(f(1))-vdW1(r,(-f(1)-1+ddelta)) )+f(6).*(f(3)-f(2)./r./r)+f(2).*f(7)./r
            f(6)
            f(7)
            -f(7)./r+f(6)./r./r-K.*f(2).*f(2)./r./r./r
            ];
    end
% -------------------------------------------------------------------------
    function res = ex5bc(f0,finf)
        res = [f0(1)+w0
            f0(2)/conarea
            f0(5)                                                                                                                                  %Non-slipping±ß½çÌõ¼ş
            f0(6)
            % f0(7)-nu*f0(6)/conarea-(1-nu)*tau
            finf(1)
            finf(2)/infinity
            finf(7)-nu*finf(6)/infinity-(1-nu)*tau
            ];
    end
    function y = ex1init(r)
        for e = 1:length(r)
            % if r(e) > r_indenter
            %    lp=sqrt(c)*1i;
            %    ln=-sqrt(c)*1i;
            %    lp=sqrt(lp);
            %    ln=sqrt(ln);
            %    ypred(e)=-2*w0/log(ln/lp)*(besselk(0,lp^0.5*r)-besselk(0,ln^0.5*r));
            % else
            ypred(e)=-w0;
            % end

        end
        y = [      ypred
            0
            0
            0
            0
            0
            0
            ];
    end
end