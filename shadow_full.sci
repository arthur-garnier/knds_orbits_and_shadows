function shadow_full(Lambda,Mass,Kerr,Newman,Image,Accretion_data)
    //This is the "full" version of the shadowing program, that doesn't use the Weierstrass function.
    //Its use is depreciated as the other version is much faster...
    cft=0.8; itermax=50; N=1200;
    if Lambda<>0 then
        c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
        Rs=2*GSI*Mass/cSI^2; J=Kerr*GSI*Mass^2/cSI; A=J/(Mass*cSI); alpha=-Accretion_data(2);
        Q=Newman*2*Mass*sqrt(%pi*e0*GSI*(1-Kerr^2)); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
        x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000;
        rs=2; rg=1; a=Kerr; T_int=0; T_ext=0; lam=0.8; chi=1+Lambda*a^2/3;
        Mrate=Accretion_data(6); Mrate=Mrate(1)*Rs*cSI^2/(Mass*2*GSI); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
        Temp=[];
        if length(Accretion_data(6))>1 then
            T_int=Accretion_data(6); T_ext=T_int(3); T_int=T_int(2);
        end

        if Kerr==0 then
            function Y=Carter_ter(tau,V)
                r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
                Dr=(1-Lambda*r^2/3)*r^2-2*r+rq;
                Drp=-4/3*Lambda*r^3+2*r-2;
                Pofr=E*r^2;
                prp=((2*E*r*Pofr-Drp*k/2)/Dr-Drp*pr^2)/r^2;
                pthp=(cos(th)*sin(th)*Lz^2/sin(th)^4)/r^2;
                rp=Dr*pr/r^2;
                thp=pth/r^2;
                php=Lz/(r^2*sin(th)^2);
                Y=[rp,thp,php,prp,pthp]';
            endfunction
        else
            function Y=Carter_ter(tau,V)
                r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
                Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                Drp=-2/3*Lambda*a^2*r-4/3*Lambda*r^3+2*r-2; Dtp=-2/3*Lambda*a^2*cos(th)*sin(th);
                Pofr=chi*(E*(r^2+a^2)-a*Lz); Wofth=chi*(a*E*sin(th)-Lz/sin(th));
                prp=((2*chi*E*r*Pofr-Drp*k/2)/Dr-Drp*pr^2)/S;
                pthp=((Dtp*k/2-chi^2*cos(th)*sin(th)*(a^2*E^2-Lz^2/sin(th)^4))/Dt-Dtp*pth^2)/S;
                rp=Dr*pr/S;
                thp=Dt*pth/S;
                php=chi/S*(a*Pofr/Dr-Wofth/(Dt*sin(th)));
                Y=[rp,thp,php,prp,pthp]';
            endfunction
        end


        function wp=projtoplane_bis(w)
            x=w(1); y=w(2); z=w(3); tt=atan(-sqrt(y^2+z^2)/x);
            wp=[-rf,y*tt*rf/sqrt(y^2+z^2),z*tt*rf/sqrt(y^2+z^2)]';
        endfunction

        function BB=BoyerLindquist_bis(R,T,P)
            BB=[sqrt(R^2+A^2)*sin(T)*cos(P),sqrt(R^2+A^2)*sin(T)*sin(P),R*cos(T)]';
        endfunction

        function XX=InvBoyerLindquist_bis(x,xp,y,yp,z,zp)
            P=atan(y,x);
            R=sqrt((-A^2+x^2+y^2+z^2+sqrt(A^2*(A^2-2*x^2-2*y^2+2*z^2)+(x^2+y^2+z^2)^2))/2);
            T=acos(z/R);
            Rp=R*(x*xp+y*yp+z*zp)/(2*R^2+A^2-x^2-y^2-z^2)+A^2*z*zp/(R*(2*R^2+A^2-x^2-y^2-z^2));
            Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
            Pp=(yp*x-xp*y)/(x^2+y^2);
            XX=[0,1,R,Rp,T,Tp,P,Pp]';
        endfunction

        Xmax=22983; tau=-2*cSI/Rs*0.00042; dtau=tau/N; 

        if Accretion_data(1)<2 then
            Img=imread(Image); IMG=[]; IMG(1,1,1)=0; Npix=size(Img); Npiy=Npix(2); Npix=Npix(1);
            for i=[1:3]
                IMG(1:Npiy,1:Npix,i)=(Img(:,:,i)');
            end
            IMG=double(IMG)/256;
            Npix=size(IMG); Npiy=Npix(2); Npix=Npix(1);
            XX=linspace(-Xmax,Xmax,Npix); YY=linspace(-Xmax*Npiy/Npix,Xmax*Npiy/Npix,Npiy);
            h=x0*Xmax*sqrt(1+Npiy^2/Npix^2)/(rf-Xmax*sqrt(1+Npiy^2/Npix^2));
        else
            Npix=Accretion_data(1); Npiy=Npix;
            XX=linspace(-Xmax,Xmax,Npix); YY=linspace(-Xmax,Xmax,Npiy);
            h=x0*Xmax*sqrt(2)/(rf-Xmax*sqrt(2));
        end

        function Z=init_conds_with_angle_bis(y,z)
            v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]'; matrot=[cos(alpha),0,-sin(alpha);0,1,0;sin(alpha),0,cos(alpha)];
            vrot=matrot*[v0(1);v0(3);v0(5)]; vvrot=-matrot*[v0(2);v0(4);v0(6)];
            Z=InvBoyerLindquist_bis(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
            Z=[Z(3);Z(5);Z(7);Z(4);Z(6);Z(8)];
            Z=[2/Rs*Z(1);Z(2);Z(3);Z(4)/cSI;Z(5)*Rs/(2*cSI);Z(6)*Rs/(2*cSI)];
        endfunction

        function vel=velocity(sphvel)
            r=sphvel(1); th=sphvel(2); ph=sphvel(3); rp=sphvel(4); thp=sphvel(5); php=sphvel(6);
            vx=cos(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) - sqrt(r^2 + a^2)*php*sin(ph)*sin(th) + sqrt(r^2 + a^2)*cos(ph)*thp*cos(th);
            vy=sin(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) + sqrt(r^2 + a^2)*php*cos(ph)*sin(th) + sqrt(r^2 + a^2)*sin(ph)*thp*cos(th);
            vz=rp*cos(th) - r*thp*sin(th);
            vel=[vx,vy,vz];
        endfunction

        function rcol=doppler_color(dope)
            dil_dope=(dope-1/2)^1;
            if dil_dope<1/4 then
                rcol=[0,4*dil_dope,1];
            elseif (1/4<=dil_dope & dil_dope<1/2) then
                rcol=[0,1,2-4*dil_dope];
            elseif (1/2<=dil_dope & dil_dope<3/4) then
                rcol=[4*dil_dope-2,1,0];
            else
                rcol=[1,4-4*dil_dope,0];
            end
            if (abs(rcol(1))<0.05 & abs(rcol(3))<0.05 & abs(rcol(2)-1)<0.05) then rcol=[2,2,2];
            end
        endfunction

        if (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/doppler_shift;
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/doppler_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/grav_shift;
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/(grav_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/(grav_shift*doppler_shift);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/(grav_shift*doppler_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)==" ") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/doppler_shift,doppler_coeff^2,1/doppler_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/grav_shift,doppler_coeff^2,1/grav_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/(grav_shift*doppler_shift)^0,doppler_coeff^2,1/(grav_shift*doppler_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)==" ") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright,doppler_coeff^2,1,T];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                colou=doppler_color(doppler_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                colou=doppler_color(grav_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(Dr-a^2*sin(th)^2*Dt)/(chi^2*S));
                colou=doppler_color(doppler_shift*grav_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        else
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); lam=0.2;
                cb=ones(1,4); cb=[255,69,0]/256; cb(4)=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
            endfunction
        end


        if Accretion_data(1)==0 then
            Xred=list(); Umax=0; Vmax=0;
            for y=XX
                for z=YY
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                    E=sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(rp^2*S/Dr+thp^2*S/Dt));
                    Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi^2*E*(r^2+a^2)-php*S*Dr)-a*chi^2*E*Dr);
                    pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
                    Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+Lambda^2/3*(a*E-Lz)^2));
                    k=Q+chi^2*(a*E-Lz)^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X];
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch')
                    if AAA==0 then Vec=ode(meth,X,0,[0:dtau:tau],Carter_ter);
                    end
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:);
                    wef=zeros(3,1);
                    for i=[1:length(R)]
                        bb=BoyerLindquist_bis(R(i),theta(i),phi(i)); bb=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*bb;
                        rb=sqrt(bb(1)^2+bb(2)^2+bb(3)^2);
                        if (abs(rb-rf)<2.5e2 & wef==zeros(3,1)) then wef=bb;
                        end
                    end
                    if wef<>zeros(3,1) then
                        wef=projtoplane_bis(wef);
                        Xred($+1)=[y,wef(2);z,wef(3)]; Umax=max(Umax,abs(wef(2))); Vmax=max(Vmax,abs(wef(3)));
                    end
                end
            end

            //Here is the supplementary loop: we computed the maximal coordinates from the rays found above
            //and using these values, we now compute the RGB value of each pixel. This slows the program by an order of magnitude.
            //Though a priori more accurate than with the Weierstrass trick, the resulting images are very similar,
            //when a relatively hign resolution is wanted (typically, more than 100 pixels).
            Umax=Umax*(1+1/1000); Vmax=Vmax*(1+1/1000); xred=[]; xred(1,1,1)=0;
            for i=[1:Npix]
                for j=[1:Npiy]
                    wef=zeros(2,2);
                    for xe=Xred
                        if (xe(1,1)==XX(i) & xe(2,1)==YY(j)) then wef=xe;
                        end
                    end
                    if wef==zeros(2,2) then xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    else s1=min(1,abs(wef(1,2)+Umax)/(2*Umax)); s2=min(1,abs(wef(2,2)+Vmax)/(2*Vmax)); ii=ceil(s1*Npix); jj=ceil(s2*Npiy); xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
            end

        elseif Accretion_data(1)==1
            Xred=list(); Umax=0; Vmax=0; dop_max=0;
            for y=XX
                for z=YY
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                    E=sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(rp^2*S/Dr+thp^2*S/Dt));
                    Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi^2*E*(r^2+a^2)-php*S*Dr)-a*chi^2*E*Dr);
                    pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
                    Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+Lambda^2/3*(a*E-Lz)^2));
                    k=Q+chi^2*(a*E-Lz)^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X]; Tau=tau;
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch','n'); iter=0;
                    while (AAA<>0 & iter<itermax)
                        Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                    end
                    Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                    vef=zeros(1,4); wef=zeros(3,1);
                    for i=[1:length(R)]
                        rb=sqrt(Vec(1,i)^2+a^2);
                        if (rb*Rs/2>rint & rb*Rs/2<rext & abs(theta(i)-%pi/2)<1/100 & vef==zeros(1,4)) then 
                            vef=[BoyerLindquist_bis(R(i),theta(i),phi(i))',accretion_disk([R(i),theta(i),phi(i),PR(i),PTH(i)])];
                            Temp=[Temp,vef($)];
                        end
                    end
                    for i=[1:length(R)]
                        bb=BoyerLindquist_bis(R(i),theta(i),phi(i));
                        bb=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*bb;
                        rbb=sqrt(bb(1)^2+bb(2)^2+bb(3)^2);
                        if (abs(rbb-rf)<2.5e2 & wef==zeros(3,1)) then wef=bb;
                        end
                    end
                    if wef<>zeros(3,1) then
                        wef=projtoplane_bis(wef);
                        Xred($+1)=[y,wef(2),1,1;z,wef(3),0,0;0,0,0,0];
                        Umax=max(Umax,abs(wef(2))); Vmax=max(Vmax,abs(wef(3)));
                    end
                    if vef<>zeros(1,4) then
                        if (Accretion_data(3)=="Custom") then
                            Xred($+1)=[y,vef(1),0,0;z,vef(2),vef(7),vef(9);vef(4),vef(5),vef(6),vef(8)];
                            dop_max=max(dop_max,vef(8));
                        else
                            Xred($+1)=[y,vef(1),0,0;z,vef(2),vef(7),0;vef(4),vef(5),vef(6),0];
                        end
                        Umax=max(Umax,abs(vef(1))); Vmax=max(Vmax,abs(vef(2)));
                    end
                end
            end
            Umax=Umax*(1+1/1000); Vmax=Vmax*(1+1/1000); xred=[]; xred(1,1,1)=0;
            for i=[1:Npix]
                for j=[1:Npiy]
                    wef=zeros(2,2);
                    for xe=Xred
                        if (xe(1,1)==XX(i) & xe(2,1)==YY(j)) then wef=xe;
                        end
                    end
                    if wef==zeros(2,2) then xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    else 
                        if wef(1,3)==1 then
                            s1=min(1,abs(wef(1,2)+Umax)/(2*Umax)); s2=min(1,abs(wef(2,2)+Vmax)/(2*Vmax)); ii=ceil(s1*Npix); jj=ceil(s2*Npiy); xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        else
                            if (Accretion_data(3)=="Custom") then
                                flo=floor(wef(2,4)*(T_int+(T_ext-T_int)*wef(3,4)/dop_max));
                                wef(3,1:3)=blackbody(find(abs(flo-blackbody(:,1))==min(abs(flo-blackbody(:,1))))(1),2:4);
                            end 
                            xred(i,j,1)=wef(2,3)*wef(3,1); xred(i,j,2)=wef(2,3)*wef(3,2); xred(i,j,3)=wef(2,3)*wef(3,3);
                        end
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
            end

        elseif Accretion_data(1)>1 then
            Xred=list(); Umax=0; Vmax=0; dop_max=zeros(Npix,Npiy); xred=[]; xred(1,1,1)=0;
            for i=1:Npix
                y=XX(i)
                for j=1:Npiy
                    z=YY(j);
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                    E=sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(rp^2*S/Dr+thp^2*S/Dt));
                    Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi^2*E*(r^2+a^2)-php*S*Dr)-a*chi^2*E*Dr);
                    pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
                    Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+Lambda^2/3*(a*E-Lz)^2));
                    k=Q+chi^2*(a*E-Lz)^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X]; Tau=tau;
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch','n'); iter=0;
                    while (AAA<>0 & iter<itermax)
                        Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                    end
                    Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                    xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    ll=find(abs(theta-%pi/2)<1/100); whitness=0;
                    for l=ll
                        rb=sqrt(Vec(1,l)^2+a^2);
                        if (rb*Rs/2>rint & rb*Rs/2<rext & whitness==0) then
                            whitness=1;
                            vef=[BoyerLindquist_bis(R(l),theta(l),phi(l))',accretion_disk([R(l),theta(l),phi(l),PR(l),PTH(l)])];
                            Temp=[Temp,vef($)];
                            if Accretion_data(3)=="Custom" then
                                xred(i,j,1)=-exp(1); xred(i,j,2)=vef(7); xred(i,j,3)=vef(9); dop_max(i,j)=vef(8);
                            else
                                xred(i,j,1)=vef(7)*vef(4); xred(i,j,2)=vef(7)*vef(5); xred(i,j,3)=vef(7)*vef(6);
                            end
                        end
                    end
                end
            end

            dp_max=max(dop_max);
            if Accretion_data(3)=="Custom" then
                for i=1:Npix
                    for j=1:Npiy
                        if xred(i,j,1)==-exp(1) then
                            flo=floor(xred(i,j,3)*(T_int+(T_ext-T_int)*dop_max(i,j)/dp_max));
                            wef=blackbody(find(abs(flo-blackbody(:,1))==min(abs(flo-blackbody(:,1))))(1),2:4);
                            xred(i,j,1)=xred(i,j,2)*wef(1); xred(i,j,3)=xred(i,j,2)*wef(3); xred(i,j,2)=xred(i,j,2)*wef(2);
                        end
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Accretion_data(1),1:Accretion_data(1),i)=xred(:,:,i)';
            end
        end

        figure()
        imshow(xredt);
        aa=gca();
        aa.isoview="on";
        if (Accretion_data(8)==1 & Accretion_data(3)<>" ") then
            if Accretion_data(3)=="Custom" then
                Text=T_ext; Tint=T_int;
            else
                Text=min(Temp); Tint=max(Temp);
            end
            h=gcf(); h.color_map=whitecolormap(1);
            n=find(abs(Text-blackbody(:,1))==min(abs(Text-blackbody(:,1))))(1);
            while blackbody(n,1)<Tint
                addcolor(blackbody(n,2:4)); n=n+1;
            end
            colorbar(Text,Tint); colbar = gce();
            colbar.type
            ylabel(colbar, "Temperature [K]");
            colbar.title.font_size = 3;
            aa = gcf(); aa.background = 1; aa=gce(); aa=aa.children(1); aa.background=0;
        end
    else
        c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
        Rs=2*GSI*Mass/cSI^2; J=Kerr*GSI*Mass^2/cSI; A=J/(Mass*cSI); alpha=-Accretion_data(2);
        Q=Newman*2*Mass*sqrt(%pi*e0*GSI*(1-Kerr^2)); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
        x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000;
        rs=2; rg=1; a=Kerr; T_int=0; T_ext=0; lam=0.8; chi=1;
        Mrate=Accretion_data(6); Mrate=Mrate(1)*Rs*cSI^2/(Mass*2*GSI); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
        Temp=[];
        if length(Accretion_data(6))>1 then
            T_int=Accretion_data(6); T_ext=T_int(3); T_int=T_int(2);
        end

        if Kerr<>0 then
            function Y=Carter_ter(tau,V)
                r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
                D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                prp=(k*(1-r)+2*r*E^2*(r^2+a^2)-2*a*E*Lz)/(S*D)-(2*pr^2*(r-1))/S;
                pthp=sin(th)*cos(th)*(Lz^2/sin(th)^4-a^2*E^2)/S;
                rp=D*pr/S;
                thp=pth/S;
                php=a*(E*(r^2+a^2)-a*Lz)/(S*D)+(Lz/sin(th)^2-a*E)/S;
                Y=[rp,thp,php,prp,pthp]';
            endfunction
        else
            function Y=Carter_ter(tau,V)
                r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
                D=r^2-2*r+rq; S=r^2;
                prp=(k*(1-r)+2*r*E^2*(r^2))/(S*D)-(2*pr^2*(r-1))/S;
                pthp=sin(th)*cos(th)*(Lz^2/sin(th)^4)/S;
                rp=D*pr/S;
                thp=pth/S;
                php=a*(E*(r^2+a^2)-a*Lz)/(S*D)+(Lz/sin(th)^2-a*E)/S;
                Y=[rp,thp,php,prp,pthp]';
            endfunction
        end


        function wp=projtoplane_bis(w)
            x=w(1); y=w(2); z=w(3); tt=atan(-sqrt(y^2+z^2)/x);
            wp=[-rf,y*tt*rf/sqrt(y^2+z^2),z*tt*rf/sqrt(y^2+z^2)]';
        endfunction

        function BB=BoyerLindquist_bis(R,T,P)
            BB=[sqrt(R^2+A^2)*sin(T)*cos(P),sqrt(R^2+A^2)*sin(T)*sin(P),R*cos(T)]';
        endfunction

        function XX=InvBoyerLindquist_bis(x,xp,y,yp,z,zp)
            P=atan(y,x);
            R=sqrt((-A^2+x^2+y^2+z^2+sqrt(A^2*(A^2-2*x^2-2*y^2+2*z^2)+(x^2+y^2+z^2)^2))/2);
            T=acos(z/R);
            Rp=R*(x*xp+y*yp+z*zp)/(2*R^2+A^2-x^2-y^2-z^2)+A^2*z*zp/(R*(2*R^2+A^2-x^2-y^2-z^2));
            Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
            Pp=(yp*x-xp*y)/(x^2+y^2);
            XX=[0,1,R,Rp,T,Tp,P,Pp]';
        endfunction

        Xmax=22983; tau=-2*cSI/Rs*0.00042; dtau=tau/N; 

        if Accretion_data(1)<2 then
            Img=imread(Image); IMG=[]; IMG(1,1,1)=0; Npix=size(Img); Npiy=Npix(2); Npix=Npix(1);
            for i=[1:3]
                IMG(1:Npiy,1:Npix,i)=(Img(:,:,i)');
            end
            IMG=double(IMG)/256;
            Npix=size(IMG); Npiy=Npix(2); Npix=Npix(1);
            XX=linspace(-Xmax,Xmax,Npix); YY=linspace(-Xmax*Npiy/Npix,Xmax*Npiy/Npix,Npiy);
            h=x0*Xmax*sqrt(1+Npiy^2/Npix^2)/(rf-Xmax*sqrt(1+Npiy^2/Npix^2));
        else
            Npix=Accretion_data(1); Npiy=Npix;
            XX=linspace(-Xmax,Xmax,Npix); YY=linspace(-Xmax,Xmax,Npiy);
            h=x0*Xmax*sqrt(2)/(rf-Xmax*sqrt(2));
        end

        function Z=init_conds_with_angle_bis(y,z)
            v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]'; matrot=[cos(alpha),0,-sin(alpha);0,1,0;sin(alpha),0,cos(alpha)];
            vrot=matrot*[v0(1);v0(3);v0(5)]; vvrot=-matrot*[v0(2);v0(4);v0(6)];
            Z=InvBoyerLindquist_bis(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
            Z=[Z(3);Z(5);Z(7);Z(4);Z(6);Z(8)];
            Z=[2/Rs*Z(1);Z(2);Z(3);Z(4)/cSI;Z(5)*Rs/(2*cSI);Z(6)*Rs/(2*cSI)];
        endfunction

        function vel=velocity(sphvel)
            r=sphvel(1); th=sphvel(2); ph=sphvel(3); rp=sphvel(4); thp=sphvel(5); php=sphvel(6);
            vx=cos(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) - sqrt(r^2 + a^2)*php*sin(ph)*sin(th) + sqrt(r^2 + a^2)*cos(ph)*thp*cos(th);
            vy=sin(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) + sqrt(r^2 + a^2)*php*cos(ph)*sin(th) + sqrt(r^2 + a^2)*sin(ph)*thp*cos(th);
            vz=rp*cos(th) - r*thp*sin(th);
            vel=[vx,vy,vz];
        endfunction

        function rcol=doppler_color(dope)
            dil_dope=(dope-1/2)^1;
            if dil_dope<1/4 then
                rcol=[0,4*dil_dope,1];
            elseif (1/4<=dil_dope & dil_dope<1/2) then
                rcol=[0,1,2-4*dil_dope];
            elseif (1/2<=dil_dope & dil_dope<3/4) then
                rcol=[4*dil_dope-2,1,0];
            else
                rcol=[1,4-4*dil_dope,0];
            end
            if (abs(rcol(1))<0.05 & abs(rcol(3))<0.05 & abs(rcol(2)-1)<0.05) then rcol=[2,2,2];
            end
        endfunction

        if (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/doppler_shift;
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/doppler_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/grav_shift;
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/(grav_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4); T=T/(grav_shift*doppler_shift);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright/(grav_shift*doppler_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Black-body" & Accretion_data(4)==" ") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                colou=blackbody(find(abs(T-blackbody(:,1))==min(abs(T-blackbody(:,1))))(1),2:4);
                cb=[colou,bright,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/doppler_shift,doppler_coeff^2,1/doppler_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/grav_shift,doppler_coeff^2,1/grav_shift,T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright/(grav_shift*doppler_shift)^0,doppler_coeff^2,1/(grav_shift*doppler_shift),T];
            endfunction
        elseif (Accretion_data(3)=="Custom" & Accretion_data(4)==" ") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2);
                doppler_coeff=1-sqrt((2*rb^2 + a*(a - 4)*rb + 2*a^2)/rb^3);
                T=(T0/(rb*Rs)^3*(1-sqrt(2*rint/(rb*Rs))))^(1/4);
                if Accretion_data(7)<>0 then
                    bright=Accretion_data(7)*4.086e-21*T^5;
                else
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                end
                cb=[0,0,0,bright,doppler_coeff^2,1,T];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                colou=doppler_color(doppler_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Gravitation") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                colou=doppler_color(grav_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler+") then
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
                veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));
                velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                grav_shift=1/sqrt(abs(D-a^2*sin(th)^2)/S);
                colou=doppler_color(doppler_shift*grav_shift);
                bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                cb=[colou,bright];
            endfunction
        else
            function cb=accretion_disk(V)
                r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); lam=0.2;
                cb=ones(1,4); cb=[255,69,0]/256; cb(4)=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
            endfunction
        end


        if Accretion_data(1)==0 then
            Xred=list(); Umax=0; Vmax=0;
            for y=XX
                for z=YY
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
                    E=sqrt((D-a^2*sin(th)^2)*(rp^2+D*thp^2)/D+D*php^2*sin(th)^2);
                    Lz=((S*D*php+a*E*(rq-2*r))*sin(th)^2)/(D-a^2*sin(th)^2);
                    pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
                    Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*E^2);
                    k=Q+Lz^2+a^2*E^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X];
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch')
                    if AAA==0 then Vec=ode(meth,X,0,[0:dtau:tau],Carter_ter);
                    end
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:);
                    wef=zeros(3,1);
                    for i=[1:length(R)]
                        bb=BoyerLindquist_bis(R(i),theta(i),phi(i)); bb=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*bb;
                        rb=sqrt(bb(1)^2+bb(2)^2+bb(3)^2);
                        if (abs(rb-rf)<2.5e2 & wef==zeros(3,1)) then wef=bb;
                        end
                    end
                    if wef<>zeros(3,1) then
                        wef=projtoplane_bis(wef);
                        Xred($+1)=[y,wef(2);z,wef(3)]; Umax=max(Umax,abs(wef(2))); Vmax=max(Vmax,abs(wef(3)));
                    end
                end
            end
            Umax=Umax*(1+1/1000); Vmax=Vmax*(1+1/1000); xred=[]; xred(1,1,1)=0;
            for i=[1:Npix]
                for j=[1:Npiy]
                    wef=zeros(2,2);
                    for xe=Xred
                        if (xe(1,1)==XX(i) & xe(2,1)==YY(j)) then wef=xe;
                        end
                    end
                    if wef==zeros(2,2) then xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    else s1=min(1,abs(wef(1,2)+Umax)/(2*Umax)); s2=min(1,abs(wef(2,2)+Vmax)/(2*Vmax)); ii=ceil(s1*Npix); jj=ceil(s2*Npiy); xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
            end

        elseif Accretion_data(1)==1
            Xred=list(); Umax=0; Vmax=0; dop_max=0;
            for y=XX
                for z=YY
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
                    E=sqrt((D-a^2*sin(th)^2)*(rp^2+D*thp^2)/D+D*php^2*sin(th)^2);
                    Lz=((S*D*php+a*E*(rq-2*r))*sin(th)^2)/(D-a^2*sin(th)^2);
                    pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
                    Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*E^2);
                    k=Q+Lz^2+a^2*E^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X]; Tau=tau;
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch','n'); iter=0;
                    while (AAA<>0 & iter<itermax)
                        Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                    end
                    Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                    vef=zeros(1,4); wef=zeros(3,1);
                    for i=[1:length(R)]
                        rb=sqrt(Vec(1,i)^2+a^2);
                        if (rb*Rs/2>rint & rb*Rs/2<rext & abs(theta(i)-%pi/2)<1/100 & vef==zeros(1,4)) then 
                            vef=[BoyerLindquist_bis(R(i),theta(i),phi(i))',accretion_disk([R(i),theta(i),phi(i),PR(i),PTH(i)])];
                            Temp=[Temp,vef($)];
                        end
                    end
                    for i=[1:length(R)]
                        bb=BoyerLindquist_bis(R(i),theta(i),phi(i));
                        bb=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*bb;
                        rbb=sqrt(bb(1)^2+bb(2)^2+bb(3)^2);
                        if (abs(rbb-rf)<2.5e2 & wef==zeros(3,1)) then wef=bb;
                        end
                    end
                    if wef<>zeros(3,1) then
                        wef=projtoplane_bis(wef);
                        Xred($+1)=[y,wef(2),1,1;z,wef(3),0,0;0,0,0,0];
                        Umax=max(Umax,abs(wef(2))); Vmax=max(Vmax,abs(wef(3)));
                    end
                    if vef<>zeros(1,4) then
                        if (Accretion_data(3)=="Custom") then
                            Xred($+1)=[y,vef(1),0,0;z,vef(2),vef(7),vef(9);vef(4),vef(5),vef(6),vef(8)];
                            dop_max=max(dop_max,vef(8));
                        else
                            Xred($+1)=[y,vef(1),0,0;z,vef(2),vef(7),0;vef(4),vef(5),vef(6),0];
                        end
                        Umax=max(Umax,abs(vef(1))); Vmax=max(Vmax,abs(vef(2)));
                    end
                end
            end
            Umax=Umax*(1+1/1000); Vmax=Vmax*(1+1/1000); xred=[]; xred(1,1,1)=0;
            for i=[1:Npix]
                for j=[1:Npiy]
                    wef=zeros(2,2);
                    for xe=Xred
                        if (xe(1,1)==XX(i) & xe(2,1)==YY(j)) then wef=xe;
                        end
                    end
                    if wef==zeros(2,2) then xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    else 
                        if wef(1,3)==1 then
                            s1=min(1,abs(wef(1,2)+Umax)/(2*Umax)); s2=min(1,abs(wef(2,2)+Vmax)/(2*Vmax)); ii=ceil(s1*Npix); jj=ceil(s2*Npiy); xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        else
                            if (Accretion_data(3)=="Custom") then
                                flo=floor(wef(2,4)*(T_int+(T_ext-T_int)*wef(3,4)/dop_max));
                                wef(3,1:3)=blackbody(find(abs(flo-blackbody(:,1))==min(abs(flo-blackbody(:,1))))(1),2:4);
                            end 
                            xred(i,j,1)=wef(2,3)*wef(3,1); xred(i,j,2)=wef(2,3)*wef(3,2); xred(i,j,3)=wef(2,3)*wef(3,3);
                        end
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
            end

        elseif Accretion_data(1)>1 then
            Xred=list(); Umax=0; Vmax=0; dop_max=zeros(Npix,Npiy); xred=[]; xred(1,1,1)=0;
            for i=1:Npix
                y=XX(i)
                for j=1:Npiy
                    z=YY(j);
                    X=init_conds_with_angle_bis(y,z);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
                    E=sqrt((D-a^2*sin(th)^2)*(rp^2+D*thp^2)/D+D*php^2*sin(th)^2);
                    Lz=((S*D*php+a*E*(rq-2*r))*sin(th)^2)/(D-a^2*sin(th)^2);
                    pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
                    Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*E^2);
                    k=Q+Lz^2+a^2*E^2;
                    X=[r,th,ph,pr,pth]';
                    Vec=[X]; Tau=tau;
                    AAA=execstr('ode(meth,X,0,[0:dtau:tau],Carter_ter)','errcatch','n'); iter=0;
                    while (AAA<>0 & iter<itermax)
                        Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                    end
                    Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                    R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                    xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                    ll=find(abs(theta-%pi/2)<1/100); whitness=0;
                    for l=ll
                        rb=sqrt(Vec(1,l)^2+a^2);
                        if (rb*Rs/2>rint & rb*Rs/2<rext & whitness==0) then
                            whitness=1;
                            vef=[BoyerLindquist_bis(R(l),theta(l),phi(l))',accretion_disk([R(l),theta(l),phi(l),PR(l),PTH(l)])];
                            Temp=[Temp,vef($)];
                            if Accretion_data(3)=="Custom" then
                                xred(i,j,1)=-exp(1); xred(i,j,2)=vef(7); xred(i,j,3)=vef(9); dop_max(i,j)=vef(8);
                            else
                                xred(i,j,1)=vef(7)*vef(4); xred(i,j,2)=vef(7)*vef(5); xred(i,j,3)=vef(7)*vef(6);
                            end
                        end
                    end
                end
            end
            dp_max=max(dop_max);
            if Accretion_data(3)=="Custom" then
                for i=1:Npix
                    for j=1:Npiy
                        if xred(i,j,1)==-exp(1) then
                            flo=floor(xred(i,j,3)*(T_int+(T_ext-T_int)*dop_max(i,j)/dp_max));
                            wef=blackbody(find(abs(flo-blackbody(:,1))==min(abs(flo-blackbody(:,1))))(1),2:4);
                            xred(i,j,1)=xred(i,j,2)*wef(1); xred(i,j,3)=xred(i,j,2)*wef(3); xred(i,j,2)=xred(i,j,2)*wef(2);
                        end
                    end
                end
            end
            xredt=[]; xredt(1,1,1)=0;
            for i=[1:3]
                xredt(1:Accretion_data(1),1:Accretion_data(1),i)=xred(:,:,i)';
            end
        end

        figure()
        imshow(xredt);
        aa=gca();
        aa.isoview="on";
        if (Accretion_data(8)==1 & Accretion_data(3)<>" ") then
            if Accretion_data(3)=="Custom" then
                Text=T_ext; Tint=T_int;
            else
                Text=min(Temp); Tint=max(Temp);
            end
            h=gcf(); h.color_map=whitecolormap(1);
            n=find(abs(Text-blackbody(:,1))==min(abs(Text-blackbody(:,1))))(1);
            while blackbody(n,1)<Tint
                addcolor(blackbody(n,2:4)); n=n+1;
            end
            colorbar(Text,Tint); colbar = gce();
            colbar.type
            ylabel(colbar, "Temperature [K]");
            colbar.title.font_size = 3;
            aa = gcf(); aa.background = 1; aa=gce(); aa=aa.children(1); aa.background=0;
        end
    end
endfunction
