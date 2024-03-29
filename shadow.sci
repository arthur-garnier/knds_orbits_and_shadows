//The shadowing program is divided in two subprograms:
//the first one solves the Carter equations while the second one takes advantage of the case where a=0 and uses Weirstrass' elliptic functions instead:
function xredt=shadow(Lambda,Mass,Kerr,Newman,Image,Accretion_data)
    if Kerr<>0 then //Carter's equations
        //We define the refining parameters for details near the horizon(s)
        //cft=0.95; itermax=100; N=1500;//N=3000;
        cft=0.85; itermax=50; N=1500;//N=1200;
        //cft=0.8; itermax=50; N=1200;
        //Here again, distinguish between whether Lambda is zero or not (faster if Lambda=0)    
        if Lambda<>0 then
            //Initialize: unitless initial conditions and parameters (of the black hole and the accretion disk: size and extremal temperatures if input)
            c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
            Rs=2*GSI*Mass/cSI^2; A=Kerr*Rs/2; alpha=-Accretion_data(2);//J=Kerr*GSI*Mass^2/cSI; A=J/(Mass*cSI); alpha=-Accretion_data(2);
            rq=Newman^2;//Q=Newman*2*Mass*sqrt(%pi*e0*GSI); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
            x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000;
            rs=2; rg=1; a=Kerr; T_int=0; T_ext=0; lam=0.8; chi=1+Lambda*a^2/3;
            Mrate=Accretion_data(6); Mrate=Mrate(1)*Rs*cSI^2/(Mass*2*GSI); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
            Temp=[];
            if length(Accretion_data(6))>1 then
                T_int=Accretion_data(6); T_ext=T_int(3); T_int=T_int(2);
            end

            //Defining the 'ode' functions (distinguishing between rotating and non-rotating black hole)
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

            //Equi-rectangular projection from sphere to tangent plane
            function wp=projtoplane_bis(w)
                wp=[-1,atan(w(2),-w(1)),%pi/2-acos(w(3)/rf)]';
            endfunction

            //From and to BL coordinates (with velocities)
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

            //Accretion_data(1)<2 means that we want the image with accretion disk and otherwise, we just want the accretion disk and no image is required.
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

            //The initial datum from a pixel (y,z) on the screen (a ray leaving it with the appropriate velocity)
            function Z=init_conds_with_angle_bis(y,z)
                v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]'; matrot=[cos(alpha),0,-sin(alpha);0,1,0;sin(alpha),0,cos(alpha)];
                vrot=matrot*[v0(1);v0(3);v0(5)]; vvrot=-matrot*[v0(2);v0(4);v0(6)];
                Z=InvBoyerLindquist_bis(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
                Z=[Z(3);Z(5);Z(7);Z(4);Z(6);Z(8)];
                Z=[2/Rs*Z(1);Z(2);Z(3);Z(4)/cSI;Z(5)*Rs/(2*cSI);Z(6)*Rs/(2*cSI)];
            endfunction

            function Z=init_conds_bis(y,z)
                v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]';
                vrot=[v0(1);v0(3);v0(5)]; vvrot=-[v0(2);v0(4);v0(6)];
                Z=InvBoyerLindquist_bis(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
                Z=[Z(3);Z(5);Z(7);Z(4);Z(6);Z(8)];
                Z=[2/Rs*Z(1);Z(2);Z(3);Z(4)/cSI;Z(5)*Rs/(2*cSI);Z(6)*Rs/(2*cSI)];
            endfunction

            //Cartesian velocity from BL velocity
            function vel=velocity(sphvel)
                r=sphvel(1); th=sphvel(2); ph=sphvel(3); rp=sphvel(4); thp=sphvel(5); php=sphvel(6);
                vx=cos(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) - sqrt(r^2 + a^2)*php*sin(ph)*sin(th) + sqrt(r^2 + a^2)*cos(ph)*thp*cos(th);
                vy=sin(ph)*sin(th)*r*rp/sqrt(r^2 + a^2) + sqrt(r^2 + a^2)*php*cos(ph)*sin(th) + sqrt(r^2 + a^2)*sin(ph)*thp*cos(th);
                vz=rp*cos(th) - r*thp*sin(th);
                vel=[vx,vy,vz];
            endfunction

            //Function for the color shift
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

            //Depending on what kind of accretion is required, we create a function 'accretion_disk'
            //which computes the color to be attributed to a pixel, taking into account the various effects.
            //This function is defined at this point so that we don't need to put a selection process for each pixel:
            //the color function is chosen once and for all.
            if (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=sqrt(r^2+a^2); Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; S=r^2+a^2*cos(th)^2; Dt=1+Lambda*a^2*cos(th)^2/3;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,a,Lambda)*[pt;V(4);V(5);pph]; al=(a+r^2/sqrt(-Lambda*r^4/3+r-rq))/rb;
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
            Umax=%pi/2; Vmax=Umax*Npiy/Npix; xred=zeros(Npix,Npiy,3);

            //First case where no accretion disk is required: we compute the ray for each pixel (using Carter's equations), 
            //find if it hits the sphere and compute the associated pixel on the tangent plane.
            //Using the maximal reach previously calculated, we find its position on the projected image and give it the right color.
            //When calling the 'ode' function, we eliminate the pixels for which the integration fails (those who die in the black hole).
            if Accretion_data(1)==0 then
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
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
                        tau0=tau;//abs(1/(1+333*(Lambda)/3))*tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2))
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch');
                        if AAA==0 then Vec=ode(meth,X,0,[0:dtau:tau0],Carter_ter);
                        end
                        R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:);
                        wef=zeros(3,1); WEF=[];
                        dWEF=abs(rf-sqrt(R.^2+A^2*sin(theta).^2)); dwef=min(dWEF);
                        if dwef<2.5e2 then
                            l=find(dWEF==dwef)(1);
                            wef=BoyerLindquist_bis(R(l),theta(l),phi(l)); wef=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*wef;
                            wef=projtoplane_bis(wef);
                            s1=real(wef(2)+Umax)/(2*Umax); s2=real(wef(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
                    end
                end

                xredt=[]; xredt(1,1,1)=0;
                for i=[1:3]
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

                //If Accretion_data(1)=1, it means we want to shadow the black hole with a picture and an accretion disk.
                //We do the same as before, except that we shouldn't eliminate all the dying pixels:
                //instead, if a ray goes to the horizon (an error is then returned by 'ode'),
                //then we integrate it on [0,cft*Tau] instead of [0,Tau], where 0<cft<1 and we do it at most itermax times.
                //These parameters can be tuned at the very first line of the present function. Morally, the smaller the inner radius, 
                //the bigger cft, itermax and step size N should be.
                //The first ray that doesn't hit the horizon (if any) is kept, and we test if it crosses the theta=pi/2 plane (with a tolerance of 1e-2).
                //If so, we check if the crossing happens between the chosen radii rint and rext and if so, this adds a pixel to the accretion disk.
                //We compute the function accretion_disk for this pixel and obtain the corresponding RGB value.
            elseif Accretion_data(1)==1
                dop_max=zeros(Npix,Npiy);
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
                        X=init_conds_with_angle_bis(y,z);
                        r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                        E=sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(rp^2*S/Dr+thp^2*S/Dt));
                        Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi^2*E*(r^2+a^2)-php*S*Dr)-a*chi^2*E*Dr);
                        pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
                        Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+Lambda^2/3*(a*E-Lz)^2));
                        k=Q+chi^2*(a*E-Lz)^2;
                        X=[r,th,ph,pr,pth]';
                        tau0=tau;//abs(1/(1+333*(Lambda)/3))*tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2));
                        Vec=[X]; Tau=tau0;
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch','n'); iter=0;
                        while (AAA<>0 & iter<itermax)
                            Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                        end
                        Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                        R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                        wef=zeros(3,1); WEF=[];
                        dWEF=abs(rf-sqrt(R.^2+A^2*sin(theta).^2)); dwef=min(dWEF);
                        if dwef<2.5e2 then
                            l=find(dWEF==dwef)(1);
                            wef=BoyerLindquist_bis(R(l),theta(l),phi(l)); wef=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*wef;
                            wef=projtoplane_bis(wef);
                            s1=real(wef(2)+Umax)/(2*Umax); s2=real(wef(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
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

                //In case the inner and outer temperatures are specified (the "Custom" case), we need one more loop on pixels to find the colors of the accretion disk.
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

                //If Accretion_data(1)>1, then we don't want the picture, only the shadow of the accretion disk.
            elseif Accretion_data(1)>1 then
                dop_max=zeros(Npix,Npiy);
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
                        X=init_conds_with_angle_bis(y,z);
                        r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
                        E=sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(rp^2*S/Dr+thp^2*S/Dt));
                        Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi^2*E*(r^2+a^2)-php*S*Dr)-a*chi^2*E*Dr);
                        pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
                        Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+Lambda^2/3*(a*E-Lz)^2));
                        k=Q+chi^2*(a*E-Lz)^2;
                        X=[r,th,ph,pr,pth]';
                        tau0=tau;//abs(1/(1+333*(Lambda)/3))*tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2));
                        Vec=[X]; Tau=tau0;
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch','n'); iter=0;
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end
            end

            //Finally, we display the resulting image:
            figure()
            imshow(xredt);
            aa=gca();
            aa.isoview="on";
            if (Accretion_data(1)>0 & Accretion_data(8)==1 & Accretion_data(3)<>" ") then
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
            //The following procedure is exactly the same as above, simplified in the case where Lambda=0:
            c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
            Rs=2*GSI*Mass/cSI^2; A=Kerr*Rs/2; alpha=-Accretion_data(2);//J=Kerr*GSI*Mass^2/cSI; A=J/(Mass*cSI); alpha=-Accretion_data(2);
            rq=Newman^2;//Q=Newman*2*Mass*sqrt(%pi*e0*GSI*(1-0*Kerr^2)); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
            x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000;
            rs=2; rg=1; a=Kerr; T_int=0; T_ext=0; lam=0.8; chi=1;
            Mrate=Accretion_data(6); Mrate=Mrate(1); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
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
                wp=[-1,atan(w(2),-w(1)),%pi/2-acos(w(3)/rf)]';
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

            function Z=init_conds_bis(y,z)
                v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]';
                vrot=[v0(1);v0(3);v0(5)]; vvrot=-[v0(2);v0(4);v0(6)];
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
                dil_dope=(dope-1/2)^(1);
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
                    veloc=inv_met_mat([0,r,th,ph],2,rq,a)*[pt;V(4);V(5);pph]; al=(r^2+a*sqrt(r-rq))/sqrt((r^2+a^2)*(r-rq));;
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
            Umax=%pi/2; Vmax=Umax*Npiy/Npix; xred=zeros(Npix,Npiy,3);


            if Accretion_data(1)==0 then
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
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
                        tau0=tau;//tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2));
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch')
                        if AAA==0 then Vec=ode(meth,X,0,[0:dtau:tau0],Carter_ter);
                        end
                        R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:);
                        wef=zeros(3,1); WEF=[];
                        dWEF=abs(rf-sqrt(R.^2+A^2*sin(theta).^2)); dwef=min(dWEF);
                        if dwef<2.5e2 then
                            l=find(dWEF==dwef)(1);
                            wef=BoyerLindquist_bis(R(l),theta(l),phi(l)); wef=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*wef;
                            wef=projtoplane_bis(wef);
                            s1=real(wef(2)+Umax)/(2*Umax); s2=real(wef(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
                    end
                end

                xredt=[]; xredt(1,1,1)=0;
                for i=[1:3]
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

            elseif Accretion_data(1)==1
                dop_max=zeros(Npix,Npiy);
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
                        X=init_conds_with_angle_bis(y,z);
                        r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
                        E=sqrt((D-a^2*sin(th)^2)*(rp^2+D*thp^2)/D+D*php^2*sin(th)^2);
                        Lz=((S*D*php+a*E*(rq-2*r))*sin(th)^2)/(D-a^2*sin(th)^2);
                        pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
                        Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*E^2);
                        k=Q+Lz^2+a^2*E^2;
                        X=[r,th,ph,pr,pth]';
                        tau0=tau;//tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2));
                        Vec=[X]; Tau=tau0;
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch','n'); iter=0;
                        while (AAA<>0 & iter<itermax)
                            Tau=cft*Tau; AAA=execstr('ode(meth,X,0,[0:dtau:Tau],Carter_ter)','errcatch','n'); iter=iter+1;
                        end
                        Vec=ode(meth,X,0,[0:dtau:Tau],Carter_ter);
                        R=Rs/2*Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); PR=Vec(4,:); PTH=Vec(5,:);
                        wef=zeros(3,1); WEF=[];
                        dWEF=abs(rf-sqrt(R.^2+A^2*sin(theta).^2)); dwef=min(dWEF);
                        if dwef<2.5e2 then
                            l=find(dWEF==dwef)(1);
                            wef=BoyerLindquist_bis(R(l),theta(l),phi(l)); wef=[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)]*wef;
                            wef=projtoplane_bis(wef);
                            s1=real(wef(2)+Umax)/(2*Umax); s2=real(wef(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

            elseif Accretion_data(1)>1 then
                dop_max=zeros(Npix,Npiy);
                for y=XX
                    i=find(XX==y);
                    for z=YY
                        j=find(YY==z);
                        X=init_conds_with_angle_bis(y,z);
                        r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
                        E=sqrt((D-a^2*sin(th)^2)*(rp^2+D*thp^2)/D+D*php^2*sin(th)^2);
                        Lz=((S*D*php+a*E*(rq-2*r))*sin(th)^2)/(D-a^2*sin(th)^2);
                        pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
                        Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*E^2);
                        k=Q+Lz^2+a^2*E^2;
                        X=[r,th,ph,pr,pth]';
                        tau0=tau;//tau*(1.75-sqrt(y^2+z^2)/sqrt(max(XX)^2+max(YY)^2));
                        Vec=[X]; Tau=tau0;
                        AAA=execstr('ode(meth,X,0,[0:dtau:tau0],Carter_ter)','errcatch','n'); iter=0;
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end
            end

            figure()
            imshow(xredt);
            aa=gca();
            aa.isoview="on";
            if (Accretion_data(1)>0 & Accretion_data(8)==1 & Accretion_data(3)<>" ") then
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







        
    else //Use Weierstrass' analytic method if a=0:
        //In the case of a non-rotating black hole (a=0), the spherical symmetry of the space-time and the polar formulation of
        //the geodesic equation has an analytic solution: the Weierstrass elliptic function \wp.
        //This allows fast and accurate numerical methods (such as Newton, with suitable initialization) to find if the rays hits
        //the celestial sphere (or the accretion disk).
        if Lambda<>0 then
            //Initialize constants and some accretion data
            c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
            Rs=2*GSI*Mass/cSI^2; alpha=-Accretion_data(2); xi=Accretion_data(2); txi=tan(xi);
            x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000; rint_n=2*sizee(1); rext_n=2*sizee(2);
            rq=Newman^2;//Q=Newman*2*Mass*sqrt(%pi*e0*GSI); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
            rs=2; rg=1; a=0; T_int=0; T_ext=0; lam=0.8;
            Mrate=Accretion_data(6); Mrate=Mrate(1)*Rs*cSI^2/(Mass*2*GSI); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
            Temp=[];

            if length(Accretion_data(6))>1 then
                T_int=Accretion_data(6); T_ext=T_int(3); T_int=T_int(2);
            end

            //The Carlson algorithm for computing elliptic integrals
            function app=carlson(x,y,z)
                rtol=1e-10; xn=x; yn=y; zn=z; A0=(xn+yn+zn)/3; m=0;
                Q=46.42*max(abs(A0-xn),abs(A0-yn),abs(A0-zn)); A=A0;//Q=1/nthroot(3*rtol,6)*...
                while Q/4^m>abs(A)
                    sqx=sqrt(xn); sqy=sqrt(yn); sqz=sqrt(zn);
                    if real(sqx)<0 then
                        sqx=-sqx;
                    end
                    if real(sqy)<0 then
                        sqy=-sqy;
                    end
                    if real(sqz)<0 then
                        sqz=-sqz;
                    end
                    lm=sqx*sqy+sqx*sqz+sqy*sqz; A=(A+lm)/4;
                    xn=(xn+lm)/4; yn=(yn+lm)/4; zn=(zn+lm)/4; m=m+1;
                end
                X=(A0-x)/(4^m*A); Y=(A0-y)/(4^m*A); Z=-X-Y;
                E2=X*Y-Z^2; E3=X*Y*Z;
                app=(1-E2/10+E3/14+E2^2/24-3*E2*E3/44)/sqrt(A);
            endfunction

            //From and to spherical coordinates (with velocities), as well as the rotation function
            function BB=From_spherical(R,Rp,T,Tp,P,Pp)
                x=R*sin(T)*cos(P);
                y=R*sin(T)*sin(P);
                z=R*cos(T);
                xp=(Tp*cos(T)*cos(P)*R-sin(T)*(Pp*sin(P)*R-Rp*cos(P)));
                yp=(Tp*cos(T)*sin(P)*R+sin(T)*(Pp*cos(P)*R+Rp*sin(P)));
                zp=Rp*cos(T)-R*Tp*sin(T);
                BB=[x,xp,y,yp,z,zp]';
            endfunction
            function XX=To_spherical(x,xp,y,yp,z,zp)
                P=atan(y,x);
                R=sqrt(x^2+y^2+z^2);
                T=acos(z/R);
                Rp=(x*xp+y*yp+z*zp)/R;
                Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
                Pp=(yp*x-xp*y)/(x^2+y^2);
                XX=[R,Rp,T,Tp,P,Pp]';
            endfunction
            function v=rot(axe,theta,u)
                KK=[0,-axe(3),axe(2);axe(3),0,-axe(1);-axe(2),axe(1),0]; KK=KK/norm(axe,2);
                RR=eye(3,3)+sin(theta)*KK+(1-cos(theta))*KK^2;
                v=RR*u;
            endfunction

            //Projection on the tangent plane (equi-rectangular projection)
            function wp=projtoplane_bis(w)
                wp=[-1,atan(w(2),-w(1)),%pi/2-acos(w(3)/rf)]';
            endfunction

            Xmax=22983;

            //Accretion_data(1)<2 means that we want the image with accretion disk and otherwise, we just want the accretion disk and no image is required.
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

            //The initial datum from a pixel (y,z) on the screen (a ray leaving it with the appropriate velocity)
            function Z=init_conds_bis(y,z,alph)
                v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]'; matrot=[cos(alph),0,-sin(alph);0,1,0;sin(alph),0,cos(alph)];
                vrot=matrot*[v0(1);v0(3);v0(5)]; vvrot=-matrot*[v0(2);v0(4);v0(6)];
                Z=To_spherical(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
                Z=[Z(1);Z(3);Z(5);Z(2);Z(4);Z(6)];
                Z=[2/Rs*Z(1);Z(2);Z(3);Z(4)/cSI;Z(5)*Rs/(2*cSI);Z(6)*Rs/(2*cSI)];
            endfunction

            //Cartesian velocity from BL velocity
            function vel=velocity(sphvel)
                r=sphvel(1); th=sphvel(2); ph=sphvel(3); rp=sphvel(4); thp=sphvel(5); php=sphvel(6);
                vx=cos(ph)*sin(th)*rp - r*php*sin(ph)*sin(th) + r*cos(ph)*thp*cos(th);
                vy=sin(ph)*sin(th)*rp + r*php*cos(ph)*sin(th) + r*sin(ph)*thp*cos(th);
                vz=rp*cos(th) - r*thp*sin(th);
                vel=[vx,vy,vz];
            endfunction

            //Function for the color shift
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

            //Depending on what kind of accretion is required, we create a function 'accretion_disk'
            //which computes the color to be attributed to a pixel, taking into account the various effects.
            //This function is defined at this point so that we don't need to put a selection process for each pixel:
            //the color function is chosen once and for all.
            if (Accretion_data(3)=="Black-body" & Accretion_data(4)=="Doppler") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    grav_shift=r/sqrt(abs(DDr));
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=r/sqrt(abs(DDr));
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    grav_shift=r/sqrt(abs(DDr));
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=r/sqrt(abs(DDr));
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r;
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(abs(1-1/al^2));
                    colou=doppler_color(doppler_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            elseif (Accretion_data(3)==" " & Accretion_data(4)=="Gravitation") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    grav_shift=r/sqrt(abs(DDr));
                    colou=doppler_color(grav_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler+") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; DDr=r^2*(1-Lambda*r^2/3)-2*r+rq;
                    veloc=cosmo_inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=sqrt(r/(1-Lambda*r^3/3-rq/r));
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=r/sqrt(abs(DDr));
                    colou=doppler_color(doppler_shift*grav_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            else
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; lam=0.2;
                    cb=ones(1,4); cb=[255,69,0]/256; cb(4)=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                endfunction
            end

            //The CGL algorithm for approaching \wp:
            function zz=weierP(g2,g3,z)
                N0=12;
                zz0=z/(2^N0); zz=1./zz0^2+g2/20*zz0.^2+g3/28*zz0.^4;
                for j=1:N0
                    zz=-2*zz+(6*zz.^2-g2/2)^2./(4*(4*zz.^3-g2*zz-g3));
                end
            endfunction

            //The Newton method to see if a ray hits the celestial sphere:
            function [sol,iter,res]=newton(g2,g3,Z,t)
                function ss=toanihil(s)
                    ss=(2*rf/Rs-rbar)*(4*real(weierP(g2,g3,Z+s))-bet/3)-alp;
                endfunction
                sgn=sign(toanihil(t)); sol=t; step=-0.02;
                while sgn*toanihil(sol)>0
                    sol=sol+step;
                end
                epss=1e-12; itermax=100; iter=0;
                while (abs(toanihil(sol))>epss & iter<itermax)
                    sol=sol-toanihil(sol)/numderivative(toanihil,sol); iter=iter+1;
                end
                res=toanihil(sol);
            endfunction

            //If no accretion disk is required:
            if Accretion_data(1)==0 then
                //For each pixel in the upper right corner (i.e. y>0, z>0, this is enough by spherical symmetry), we compute its 'theta=pi/2' version (with rotations)
                //and we find the corresponding pixel on the sphere, which we project on the plane and record its coordinates.
                Xred=[];
                AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2*(1-Lambda*r^2/3)-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    //Find the constants rbar,alpha,beta,gamma,delta,g2,g3 associated to the ray:
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2+Lambda/3],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol)))(1));
                    del=(E^2)/L^2+Lambda/3; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;P];
                    end
                end

                //For each pixel in the upper-right corner, we see if it ends on the celestial sphere (using the previously computed data)
                //then, we put the position of this pixel in the list into a matrix K, which encodes the correspondence
                //between a pixel on the screen and its final state on the sphere. Besides, we extract the maximal range of coordinates.
                Umax=0; Vmax=0; Wmax=0; KK=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end
                Umax=%pi/2; Vmax=Umax*Npiy/Npix; xred=zeros(Npix,Npiy,3);

                //Attribute RGB value to the pixels
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j); xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                        if KK(i,j)<>0 then
                            P=Xred(2,KK(i,j));
                            Z=rot([1,0,0],atan(yy,xx),[cos(P);sin(P);0]);
                            Z=projtoplane_bis(rf*Z);
                            s1=(Z(2)+Umax)/(2*Umax); s2=(Z(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
                    end
                end


                xredt=[]; xredt(1,1,1)=0;
                for i=[1:3]
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

                //If accretion and picture are required:
            elseif Accretion_data(1)==1 then
                Xred=[]; Yred=[];
                AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2*(1-Lambda*r^2/3)-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2+Lambda/3],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol))));
                    del=(E^2)/L^2+Lambda/3; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;P;g2;g3;Z0;ph;alp;bet;rbar];
                    else
                        Yred(:,$+1)=[zz;0;g2;g3;Z0;ph;alp;bet;rbar];
                    end
                end

                Umax=0; Vmax=0; Wmax=0; KK=zeros(Npix,Npiy); LL=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        else
                            k=find(abs(r-Yred(1,:))==min(abs(r-Yred(1,:))))(1);
                            LL(i,j)=k; LL(i,Npiy-j+1)=k;
                            LL(Npix-i+1,j)=k; LL(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end
                Umax=%pi/2; Vmax=Umax*Npiy/Npix; xred=zeros(Npix,Npiy,3); dop_max=zeros(Npix,Npiy);

                //Now, if a pixel hits the accretion disk, we need to find the corresponding gravitational and Doppler effects.
                //To do this, we need to recover the (conjugate) momenta of the photon when it hits the disk.
                //This requires an additional loop on pixels, but the local instructions are rather easy and fast to execute and
                //this doesn't affect the execution time too much, as the whole procedure is fast enough.
                //We found convenient here to distinguish between the "upper" and "lower" parts of the disk.
                //Depending on the sign of the inclination angle, one part should be computed before the other one.
                //This is why we introduce a conditional statement on the sign of this angle.
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j); xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                        ata=atan(txi,-sign(xi)*yy/sqrt(yy^2+xx^2));
                        X=init_conds_bis(xx,yy,xi); r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        E=sqrt(rp^2+(php^2*sin(th)^2+thp^2)*(r^2*(1-Lambda*r^2/3)-2*r+rq)); L=r^2*php*sin(th)^2; pt=-E; pph=L;
                        C=r^4*thp^2+L^2*cotg(th)^2; Pth=sqrt(C);
                        if KK(i,j)==0 then
                            xe=Yred(:,LL(i,j));
                            Phi=min(%pi+ata,%pi-ata);
                            rtest=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)+(Phi-xe(6))))-xe(8)/3)+xe(9));
                            if (rint_n<rtest & rtest<rext_n & Phi<%pi/2) then
                                rtesti=rtest; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                Cf=From_spherical(Rs*rtest/2,0,%pi/2,0,Phi,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                cobra=accretion_disk([Rs*rtest/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                if Accretion_data(3)=="Custom" then
                                    xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                else
                                    xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                end
                            end
                        else
                            xe=Xred(:,KK(i,j)); 
                            Phim=%pi+atan(txi,yy/sqrt(yy^2+xx^2)); Phip=%pi-atan(txi,-yy/sqrt(yy^2+xx^2));
                            rtestp=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phip-xe(6))))-xe(8)/3)+xe(9));
                            rtestm=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phim-xe(6))))-xe(8)/3)+xe(9));
                            whit=0;
                            if alpha>=0 then
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    rtesti=rtestp; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    rtesti=rtestm; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                            else
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    rtesti=rtestm; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    rtesti=rtestp; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                            end
                            if whit==0 then
                                P=xe(2);
                                Z=rot([1,0,0],atan(yy,xx),[cos(P);sin(P);0]);
                                Z=projtoplane_bis(rf*Z);
                                s1=(Z(2)+Umax)/(2*Umax); s2=(Z(3)+Vmax)/(2*Vmax);
                                s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                                ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                                xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

                //If only the accretion disk is required:
            elseif Accretion_data(1)>1 then
                Xred=[]; Yred=[]; AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2*(1-Lambda*r^2/3)-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2+Lambda/3],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol))));
                    del=(E^2)/L^2+Lambda/3; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;2;g2;g3;Z0;ph;alp;bet;rbar];
                    else
                        Yred(:,$+1)=[zz;0;g2;g3;Z0;ph;alp;bet;rbar];
                    end
                end

                KK=zeros(Npix,Npiy); LL=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        else
                            k=find(abs(r-Yred(1,:))==min(abs(r-Yred(1,:))))(1);
                            LL(i,j)=k; LL(i,Npiy-j+1)=k;
                            LL(Npix-i+1,j)=k; LL(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end

                xred=zeros(Npix,Npiy,3); dop_max=zeros(Npix,Npiy);
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j); xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                        ata=atan(txi,-sign(xi)*yy/sqrt(yy^2+xx^2));
                        X=init_conds_bis(xx,yy,xi); r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        E=sqrt(rp^2+(php^2*sin(th)^2+thp^2)*(r^2*(1-Lambda*r^2/3)-2*r+rq)); L=r^2*php*sin(th)^2; pt=-E; pph=L;
                        C=r^4*thp^2+L^2*cotg(th)^2; Pth=sqrt(C);
                        if KK(i,j)==0 then
                            xe=Yred(:,LL(i,j));
                            Phi=min(%pi+ata,%pi-ata);
                            rtest=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)+(Phi-xe(6))))-xe(8)/3)+xe(9));
                            if (rint_n<rtest & rtest<rext_n & Phi<%pi/2) then
                                rtesti=rtest; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                Cf=From_spherical(Rs*rtest/2,0,%pi/2,0,Phi,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                cobra=accretion_disk([Rs*rtest/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                if Accretion_data(3)=="Custom" then
                                    xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                else
                                    xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                end
                            end
                        else
                            xe=Xred(:,KK(i,j)); 
                            Phim=%pi+atan(txi,yy/sqrt(yy^2+xx^2)); Phip=%pi-atan(txi,-yy/sqrt(yy^2+xx^2));
                            rtestp=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phip-xe(6))))-xe(8)/3)+xe(9));
                            rtestm=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phim-xe(6))))-xe(8)/3)+xe(9));
                            whit=0;
                            if alpha>=0 then
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    rtesti=rtestp; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    rtesti=rtestm; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                            else
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    rtesti=rtestm; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    rtesti=rtestp; DDr=rtesti^2*(1-Lambda*rtesti^2/3)-2*rtesti+rq; Pr=sqrt(E^2*rtesti^4-DDr*(L^2+C))/DDr;
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

            end

            figure()
            imshow(xredt);
            aa=gca();
            aa.isoview="on";
            if (Accretion_data(1)>0 & Accretion_data(8)==1 & Accretion_data(3)<>" ") then
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
            //The "Lambda=0" version of the above instructions:
            c=1; G=1; M=1; GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; sb=5.67e-8; meth="adams";
            Rs=2*GSI*Mass/cSI^2; alpha=-Accretion_data(2); xi=Accretion_data(2); txi=tan(xi);
            x0=50000; sizee=Accretion_data(5); rint=sizee(1)*Rs; rext=sizee(2)*Rs; rf=60000; rint_n=2*sizee(1); rext_n=2*sizee(2);
            rq=Newman^2;//Q=Newman*2*Mass*sqrt(%pi*e0*GSI); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
            rs=2; rg=1; a=0; T_int=0; T_ext=0; lam=0.8;
            Mrate=Accretion_data(6); Mrate=Mrate(1)*Rs*cSI^2/(Mass*2*GSI); T0=3*cSI^2*Rs*Mrate/(2*%pi*sb);
            Temp=[];

            if length(Accretion_data(6))>1 then
                T_int=Accretion_data(6); T_ext=T_int(3); T_int=T_int(2);
            end

            function app=carlson(x,y,z)
                rtol=1e-10; xn=x; yn=y; zn=z; A0=(xn+yn+zn)/3; m=0;
                Q=46.42*max(abs(A0-xn),abs(A0-yn),abs(A0-zn)); A=A0;//Q=1/nthroot(3*rtol,6)*...
                while Q/4^m>abs(A)
                    sqx=sqrt(xn); sqy=sqrt(yn); sqz=sqrt(zn);
                    if real(sqx)<0 then
                        sqx=-sqx;
                    end
                    if real(sqy)<0 then
                        sqy=-sqy;
                    end
                    if real(sqz)<0 then
                        sqz=-sqz;
                    end
                    lm=sqx*sqy+sqx*sqz+sqy*sqz; A=(A+lm)/4;
                    xn=(xn+lm)/4; yn=(yn+lm)/4; zn=(zn+lm)/4; m=m+1;
                end
                X=(A0-x)/(4^m*A); Y=(A0-y)/(4^m*A); Z=-X-Y;
                E2=X*Y-Z^2; E3=X*Y*Z;
                app=(1-E2/10+E3/14+E2^2/24-3*E2*E3/44)/sqrt(A);
            endfunction
            function BB=From_spherical(R,Rp,T,Tp,P,Pp)
                x=R*sin(T)*cos(P);
                y=R*sin(T)*sin(P);
                z=R*cos(T);
                xp=(Tp*cos(T)*cos(P)*R-sin(T)*(Pp*sin(P)*R-Rp*cos(P)));
                yp=(Tp*cos(T)*sin(P)*R+sin(T)*(Pp*cos(P)*R+Rp*sin(P)));
                zp=Rp*cos(T)-R*Tp*sin(T);
                BB=[x,xp,y,yp,z,zp]';
            endfunction
            function XX=To_spherical(x,xp,y,yp,z,zp)
                P=atan(y,x);
                R=sqrt(x^2+y^2+z^2);
                T=acos(z/R);
                Rp=(x*xp+y*yp+z*zp)/R;
                Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
                Pp=(yp*x-xp*y)/(x^2+y^2);
                XX=[R,Rp,T,Tp,P,Pp]';
            endfunction
            function v=rot(axe,theta,u)
                KK=[0,-axe(3),axe(2);axe(3),0,-axe(1);-axe(2),axe(1),0]; KK=KK/norm(axe,2);
                RR=eye(3,3)+sin(theta)*KK+(1-cos(theta))*KK^2;
                v=RR*u;
            endfunction

            function wp=projtoplane_bis(w)
                wp=[-1,atan(w(2),-w(1)),%pi/2-acos(w(3)/rf)]';
            endfunction

            Xmax=22983;

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

            function Z=init_conds_bis(y,z,alph)
                v0=[x0,-cSI*h/sqrt(h^2+y^2+z^2),y,cSI*y/sqrt(h^2+y^2+z^2),z,cSI*z/sqrt(h^2+y^2+z^2)]';
                matrot=[cos(alph),0,-sin(alph);0,1,0;sin(alph),0,cos(alph)];
                vrot=matrot*[v0(1);v0(3);v0(5)]; vvrot=-matrot*[v0(2);v0(4);v0(6)];
                Z=To_spherical(vrot(1),vvrot(1),vrot(2),vvrot(2),vrot(3),vvrot(3));
                Z=[Z(1);Z(3);Z(5);Z(2);Z(4);Z(6)];
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    grav_shift=1/sqrt(abs(D)/r^2);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=1/sqrt(abs(D)/r^2);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r;
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    grav_shift=1/sqrt(abs(D)/r^2);
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=1/sqrt(abs(D)/r^2);
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r;
                    doppler_coeff=1-sqrt(2/rb);
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
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(abs(1-1/al^2));
                    colou=doppler_color(doppler_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            elseif (Accretion_data(3)==" " & Accretion_data(4)=="Gravitation") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    grav_shift=1/sqrt(abs(D)/r^2);
                    colou=doppler_color(grav_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            elseif (Accretion_data(3)==" " & Accretion_data(4)=="Doppler+") then
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; D=r^2-2*r+rq;
                    veloc=inv_met_mat([0,r,th,ph],2,rq,0)*[pt;V(4);V(5);pph]; al=r/sqrt(r-rq);
                    velockep=[-sin(ph),cos(ph),0]/al; veloc=velocity([r,th,ph,veloc(2),veloc(3),veloc(4)]);
                    doppler_shift=(1-sum(veloc.*velockep)/norm(veloc,2))/sqrt(1-1/al^2);
                    grav_shift=1/sqrt(abs(D)/r^2);
                    colou=doppler_color(doppler_shift*grav_shift);
                    bright=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                    cb=[colou,bright];
                endfunction
            else
                function cb=accretion_disk(V)
                    r=2*V(1)/Rs; th=V(2); ph=V(3); rb=r; lam=0.2;
                    cb=ones(1,4); cb=[255,69,0]/256; cb(4)=1+(rb*Rs/2-rint)*(lam-1)/(rext-rint);
                endfunction
            end


            function zz=weierP(g2,g3,z)
                N0=12;
                zz0=z/(2^N0); zz=1./zz0^2+g2/20*zz0.^2+g3/28*zz0.^4;
                for j=1:N0
                    zz=-2*zz+(6*zz.^2-g2/2)^2./(4*(4*zz.^3-g2*zz-g3));
                end
            endfunction

            function [sol,iter,res]=newton(g2,g3,Z,t)
                function ss=toanihil(s)
                    ss=(2*rf/Rs-rbar)*(4*real(weierP(g2,g3,Z+s))-bet/3)-alp;
                endfunction
                sgn=sign(toanihil(t)); sol=t; step=-0.02;
                while sgn*toanihil(sol)>0
                    sol=sol+step;
                end
                epss=1e-12; itermax=100; iter=0;
                while (abs(toanihil(sol))>epss & iter<itermax)
                    sol=sol-toanihil(sol)/numderivative(toanihil,sol); iter=iter+1;
                end
                res=toanihil(sol);
            endfunction

            if Accretion_data(1)==0 then

                Xred=[];
                AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol)))(1));
                    del=(E^2)/L^2; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;P];
                    end
                end

                Umax=0; Vmax=0; Wmax=0; KK=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end
                Umax=%pi/2; Vmax=Umax*Npiy/Npix;
                xred=zeros(Npix,Npiy,3);
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j);
                        if KK(i,j)<>0 then
                            P=Xred(2,KK(i,j));
                            Z=rot([1,0,0],atan(yy,xx),[cos(P);sin(P);0]);
                            Z=projtoplane_bis(rf*Z);
                            s1=(Z(2)+Umax)/(2*Umax); s2=(Z(3)+Vmax)/(2*Vmax);
                            s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                            ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                            xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
                        end
                    end
                end


                xredt=[]; xredt(1,1,1)=0;
                for i=[1:3]
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end



            elseif Accretion_data(1)==1 then

                Xred=[]; Yred=[];
                AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol))));
                    del=(E^2)/L^2; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;P;g2;g3;Z0;ph;alp;bet;rbar];;
                    else
                        Yred(:,$+1)=[zz;0;g2;g3;Z0;ph;alp;bet;rbar];
                    end
                end

                Umax=0; Vmax=0; Wmax=0; KK=zeros(Npix,Npiy); LL=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        else
                            k=find(abs(r-Yred(1,:))==min(abs(r-Yred(1,:))))(1);
                            LL(i,j)=k; LL(i,Npiy-j+1)=k;
                            LL(Npix-i+1,j)=k; LL(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end
                Umax=%pi/2; Vmax=Umax*Npiy/Npix; xred=zeros(Npix,Npiy,3); dop_max=zeros(Npix,Npiy);
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j); xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                        ata=atan(txi,-sign(xi)*yy/sqrt(yy^2+xx^2));
                        X=init_conds_bis(xx,yy,xi); r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        E=sqrt(rp^2+(php^2*sin(th)^2+thp^2)*(r^2-2*r+rq)); L=r^2*php*sin(th)^2; pt=-E; pph=L;
                        C=r^4*thp^2+L^2*cotg(th)^2; Pth=sqrt(C);
                        if KK(i,j)==0 then
                            xe=Yred(:,LL(i,j));
                            Phi=min(%pi+ata,%pi-ata);
                            rtest=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)+(Phi-xe(6))))-xe(8)/3)+xe(9));
                            if (rint_n<rtest & rtest<rext_n & Phi<%pi/2) then
                                Pr=sqrt(rtest^2*E^2*(rtest^2+rq)-(rtest^2-2*rtest+rq)*(L^2+C))/(rtest^2-2*rtest+rq);
                                Cf=From_spherical(Rs*rtest/2,0,%pi/2,0,Phi,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                cobra=accretion_disk([Rs*rtest/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                if Accretion_data(3)=="Custom" then
                                    xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                else
                                    xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                end
                            end
                        else
                            xe=Xred(:,KK(i,j)); 
                            Phim=%pi+atan(txi,yy/sqrt(yy^2+xx^2)); Phip=%pi-atan(txi,-yy/sqrt(yy^2+xx^2));
                            rtestp=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phip-xe(6))))-xe(8)/3)+xe(9));
                            rtestm=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phim-xe(6))))-xe(8)/3)+xe(9));
                            whit=0;
                            if alpha>=0 then
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    Pr=sqrt(rtestp^2*E^2*(rtestp^2+rq)-(rtestp^2-2*rtestp+rq)*(L^2+C))/(rtestp^2-2*rtestp+rq);
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    Pr=sqrt(rtestm^2*E^2*(rtestm^2+rq)-(rtestm^2-2*rtestm+rq)*(L^2+C))/(rtestm^2-2*rtestm+rq);
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                            else
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    Pr=sqrt(rtestm^2*E^2*(rtestm^2+rq)-(rtestm^2-2*rtestm+rq)*(L^2+C))/(rtestm^2-2*rtestm+rq);
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    Pr=sqrt(rtestp^2*E^2*(rtestp^2+rq)-(rtestp^2-2*rtestp+rq)*(L^2+C))/(rtestp^2-2*rtestp+rq);
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                    whit=1;
                                end
                            end

                            if whit==0 then
                                P=xe(2);
                                Z=rot([1,0,0],atan(yy,xx),[cos(P);sin(P);0]);
                                Z=projtoplane_bis(rf*Z);
                                s1=(Z(2)+Umax)/(2*Umax); s2=(Z(3)+Vmax)/(2*Vmax);
                                s1=abs(1-abs(1-s1)); s2=abs(1-abs(1-s2));
                                ii=max(1,min(Npix,ceil(s1*Npix))); jj=max(1,min(Npiy,ceil(s2*Npiy)));
                                xred(i,j,1)=IMG(ii,jj,1); xred(i,j,2)=IMG(ii,jj,2); xred(i,j,3)=IMG(ii,jj,3);
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end



            elseif Accretion_data(1)>1 then

                Xred=[]; Yred=[]; AR=[];
                for xx=XX(floor(Npix/2)+1:$)
                    for yy=YY(floor(Npiy/2)+1:$)
                        AR=[AR,sqrt(xx^2+yy^2)];
                    end
                end
                AR=gsort(unique(AR),'c','i');

                for zz=AR
                    X=init_conds_bis(zz,0,0);
                    r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                    E=sqrt(rp^2+php^2*(r^2-2*r+rq)); L=r^2*php; pt=-E; pph=L;
                    rpol=roots(poly([-rq,2,-1,0,(E^2)/L^2],'x','c'));
                    frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
                    rbar=frpol(find(abs(frpol)==min(abs(frpol))));
                    del=(E^2)/L^2; gam=2*(2*del*rbar); bet=-1+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
                    g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
                    rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(r-rbar))+bet/12;
                    if abs(rp)<1e-12 then
                        Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    else
                        Z0=sign(-rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3));
                    end
                    [new,iter,res]=newton(g2,g3,Z0,0); P=ph+sign(php)*new;
                    if P>%pi/2 then
                        Xred(:,$+1)=[zz;2;g2;g3;Z0;ph;alp;bet;rbar];
                    else
                        Yred(:,$+1)=[zz;0;g2;g3;Z0;ph;alp;bet;rbar];
                    end
                end

                KK=zeros(Npix,Npiy); LL=zeros(Npix,Npiy);
                for i=(floor(Npix/2)+1):Npix
                    x=XX(i);
                    for j=(floor(Npiy/2)+1):Npiy
                        y=YY(j); r=sqrt(x^2+y^2);
                        if min(abs(r-Xred(1,:)))<1e-10 then
                            k=find(abs(r-Xred(1,:))==min(abs(r-Xred(1,:))))(1);
                            KK(i,j)=k; KK(i,Npiy-j+1)=k;
                            KK(Npix-i+1,j)=k; KK(Npix-i+1,Npiy-j+1)=k;
                        else
                            k=find(abs(r-Yred(1,:))==min(abs(r-Yred(1,:))))(1);
                            LL(i,j)=k; LL(i,Npiy-j+1)=k;
                            LL(Npix-i+1,j)=k; LL(Npix-i+1,Npiy-j+1)=k;
                        end
                    end
                end

                xred=zeros(Npix,Npiy,3); dop_max=zeros(Npix,Npiy);
                for i=1:Npix
                    xx=XX(i);
                    for j=1:Npiy
                        yy=YY(j); xred(i,j,1)=0; xred(i,j,2)=0; xred(i,j,3)=0;
                        ata=atan(txi,-sign(xi)*yy/sqrt(yy^2+xx^2));
                        X=init_conds_bis(xx,yy,xi); r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
                        E=sqrt(rp^2+(php^2*sin(th)^2+thp^2)*(r^2-2*r+rq)); L=r^2*php*sin(th)^2; pt=-E; pph=L;
                        C=r^4*thp^2+L^2*cotg(th)^2; Pth=sqrt(C);
                        if KK(i,j)==0 then
                            xe=Yred(:,LL(i,j));
                            Phi=min(%pi+ata,%pi-ata);
                            rtest=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)+(Phi-xe(6))))-xe(8)/3)+xe(9));
                            if (rint_n<rtest & rtest<rext_n & Phi<%pi/2) then
                                Pr=sqrt(rtest^2*E^2*(rtest^2+rq)-(rtest^2-2*rtest+rq)*(L^2+C))/(rtest^2-2*rtest+rq);
                                Cf=From_spherical(Rs*rtest/2,0,%pi/2,0,Phi,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                cobra=accretion_disk([Rs*rtest/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                if Accretion_data(3)=="Custom" then
                                    xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                else
                                    xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                end
                            end
                        else
                            xe=Xred(:,KK(i,j)); 
                            Phim=%pi+atan(txi,yy/sqrt(yy^2+xx^2)); Phip=%pi-atan(txi,-yy/sqrt(yy^2+xx^2));
                            rtestp=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phip-xe(6))))-xe(8)/3)+xe(9));
                            rtestm=(xe(7)/(4*real(weierP(xe(3),xe(4),xe(5)-(Phim-xe(6))))-xe(8)/3)+xe(9));
                            whit=0;
                            if alpha>=0 then
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    Pr=sqrt(rtestp^2*E^2*(rtestp^2+rq)-(rtestp^2-2*rtestp+rq)*(L^2+C))/(rtestp^2-2*rtestp+rq);
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    Pr=sqrt(rtestm^2*E^2*(rtestm^2+rq)-(rtestm^2-2*rtestm+rq)*(L^2+C))/(rtestm^2-2*rtestm+rq);
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                            else
                                if (rint_n<rtestm & rtestm<rext_n) then
                                    Pr=sqrt(rtestm^2*E^2*(rtestm^2+rq)-(rtestm^2-2*rtestm+rq)*(L^2+C))/(rtestm^2-2*rtestm+rq);
                                    Cf=From_spherical(Rs*rtestm/2,0,%pi/2,0,Phim,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestm/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
                                end
                                if (rint_n<rtestp & rtestp<rext_n) then
                                    Pr=sqrt(rtestp^2*E^2*(rtestp^2+rq)-(rtestp^2-2*rtestp+rq)*(L^2+C))/(rtestp^2-2*rtestp+rq);
                                    Cf=From_spherical(Rs*rtestp/2,0,%pi/2,0,Phip,0); Cf=rot([1,0,0],atan(yy,xx),[Cf(1);Cf(3);Cf(5)]);
                                    Cf=rot([0,1,0],xi,Cf); Cf=To_spherical(Cf(1),0,Cf(2),0,Cf(3),0);
                                    cobra=accretion_disk([Rs*rtestp/2,%pi/2,Cf(5),Pr,Pth]); Temp=[Temp,cobra($)];
                                    if Accretion_data(3)=="Custom" then
                                        xred(i,j,1)=-exp(1); xred(i,j,2)=cobra(4); xred(i,j,3)=cobra(6); dop_max(i,j)=cobra(5);
                                    else
                                        xred(i,j,1)=cobra(4)*cobra(1); xred(i,j,2)=cobra(4)*cobra(2); xred(i,j,3)=cobra(4)*cobra(3);
                                    end
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
                    xredt(1:Npiy,1:Npix,i)=xred(:,:,i)';
                end

            end

            figure()
            imshow(xredt);
            aa=gca();
            aa.isoview="on";
            if (Accretion_data(1)>0 & Accretion_data(8)==1 & Accretion_data(3)<>" ") then
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
    end
endfunction
