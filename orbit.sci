//The function that computes the trajectory of a test particle. 
function [Vecc,HAM,CAR]=orbit(Lambda,Mass,Kerr,Newman,IniConds,Form,Tau,N,Mu,Conserv)
    if Lambda<>0 then
        //Fundamental constants, normalized parameters and initial conditions
        GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; e=0;
        //Mass=cSI^2*Rs/(2*GSI); J=Kerr*GSI*Mass^2/cSI; a=J/(Mass*cSI); 
        Rs=2*GSI*Mass/cSI^2; J=Kerr*GSI*Mass^2/cSI; a=J/(Mass*cSI); 
        Q=Newman*2*Mass*sqrt(%pi*e0*GSI*(1-Kerr^2)); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
        G=1; c=1; Mass=1; rs=2; rg=1; a=Kerr; Q=sqrt(rq*4*%pi*e0); tau=2*cSI/Rs*Tau; chi=1+Lambda*a^2/3;
        IC=[2/Rs*IniConds(1);IniConds(2);IniConds(3);IniConds(4)/cSI;IniConds(5)*Rs/(2*cSI);IniConds(6)*Rs/(2*cSI)];
        //Compute \dot{t}_0 using Carter's formulas.
        dtau=tau/N; r0=IC(1); th0=IC(2); ph0=IC(3); rp0=IC(4); thp0=IC(5); php0=IC(6);
        Dr0=(1-Lambda*r0^2/3)*(r0^2+a^2)-2*r0+rq; Dt0=1+Lambda*a^2*cos(th0)^2/3; S0=r0^2+a^2*cos(th0)^2;
        E0=-e*r0*sqrt(rq)/(chi*S0)+sqrt(Dr0*Dt0*php0^2*sin(th0)^2/chi^4+(Dr0-Dt0*a^2*sin(th0)^2)/(chi^2*S0)*(Mu+rp0^2*S0/Dr0+thp0^2*S0/Dt0));
        L0=sin(th0)^2/(chi^2*(Dt0*a^2*sin(th0)^2-Dr0))*(Dt0*(a*chi*(chi*E0*(r0^2+a^2)+e*r0*sqrt(rq))-php0*S0*Dr0)-a*chi^2*E0*Dr0);
        time_coeff=chi/S0*((r0^2+a^2)*(chi*(E0*(r0^2+a^2)-a*L0)+e*r0*sqrt(rq))/Dr0-a*sin(th0)*(a*E0*sin(th0)-L0/sin(th0))/Dt0);
        //Select the formulation:
        if Form=="Polar" then
            Vecc=[];
            //From spherical coordinates (with velocities) to Cartesian coordinates (with velocities)
            function BB=From_spherical(R,Rp,T,Tp,P,Pp)
                x=R*sin(T)*cos(P);
                y=R*sin(T)*sin(P);
                z=R*cos(T);
                xp=(Tp*cos(T)*cos(P)*R-sin(T)*(Pp*sin(P)*R-Rp*cos(P)));
                yp=(Tp*cos(T)*sin(P)*R+sin(T)*(Pp*cos(P)*R+Rp*sin(P)));
                zp=Rp*cos(T)-R*Tp*sin(T);
                BB=[x,xp,y,yp,z,zp]';
            endfunction
            //Back to spherical coordinates (with velocities)
            function XX=To_spherical(x,xp,y,yp,z,zp)
                P=atan(y,x)
                R=sqrt(x^2+y^2+z^2);
                T=acos(z/R);
                Rp=(x*xp+y*yp+z*zp)/R;
                Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
                Pp=(yp*x-xp*y)/(x^2+y^2);
                XX=[R,Rp,T,Tp,P,Pp]';
            endfunction
            //Compute a rotated vector using Rodrigues' formula. The input is the directed axis, the angle and the vector to be rotated.
            function v=rot(axe,theta,u)
                KK=[0,-axe(3),axe(2);axe(3),0,-axe(1);-axe(2),axe(1),0]; KK=KK/norm(axe,2);
                RR=eye(3,3)+sin(theta)*KK+(1-cos(theta))*KK^2;
                v=RR*u;
            endfunction
            //Bring the initial point and velocity in the plane theta=pi/2
            X=IC; r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6); mu=Mu;
            BB=From_spherical(r,rp,th,thp,ph,php); BBi=[BB(1);BB(3);BB(5)]; BBv=[BB(2);BB(4);BB(6)];
            x=BBi(1); y=BBi(2); z=BBi(3); vx=BBv(1); vy=BBv(2); vz=BBv(3);
            if (abs(z)<1e-10 & abs(vz)<1e-10) then
                th0=0; O0=[1,0,0];
            elseif (abs(z)<1e-10 & abs(vz)>=1e-10) then
                P0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),vz]; Q0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),0];
                th0=sign(vz)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0))); O0=BBi;
            elseif (abs(z)>=1e-10 & abs(vz)<1e-10) then
                th0=%pi/2; O0=BBv;
            else
                O0=[x-z*vx/vz,y-z*vy/vz,0]; P0=[vy*z/vz-y,-vx*z/vz+x,0];
                Q0=[-z*(-vy*z+vz*y)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z*(-vx*z+vz*x)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z];
                th0=sign(z)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0)));
            end
            BBi=rot(O0,-th0,BBi); BBv=rot(O0,-th0,BBv);
            CC=To_spherical(BBi(1),BBv(1),BBi(2),BBv(2),BBi(3),BBv(3));
            R=CC(1); Rp=CC(2); T=CC(3); Tp=CC(4); P=CC(5); Pp=CC(6);
            //Angular momentum (constant) and initial polar datum
            L=R^2*Pp; X=[1/R;-Rp/L];
            //Function for 'ode'
            function Y=f(U)
                Y=[U(2);U(1)*(-2*rq*U(1)^2+3*U(1)-1-mu*rq/L^2)+mu/L^2-Lambda*mu/(3*L^2*U(1)^3)];
            endfunction
            function U=F(t,V)
                U=f(V);
            endfunction
            //Solve using 'ode'
            Vec=ode("adams",X,0,[0:dtau*Pp:tau*Pp],F);
            Vec=[1./Vec(1,:);P+Pp*[0:dtau:tau]]; Vecc=[]; HAM=[];
            test=0; ii=1;
            //Check if the trajectoy doesn't go too far or too close
            while (test==0 & ii<size(Vec)(2))
                Rf=Vec(1,ii); Pf=Vec(2,ii);
                if (Rf<0 | Rf>200*rs) then
                    test=1;
                end
                //back to the original plane and to SI units
                Cf=rot(O0,th0,[Rf*cos(Pf);Rf*sin(Pf);0]);
                Vecc=[Vecc,[Rs/2*Rf;acos(Cf(3)/Rf);atan(Cf(2),Cf(1))]]; ii=ii+1;
            end
        elseif Form=="Euler-Lagrange" then
            //initialize
            X0=[0,time_coeff,IC(1),IC(4),IC(2),IC(5),IC(3),IC(6)]'; Vec=[X0]; X=X0; CAR=[];
            //choose the right 'ode' function from parameters
            if (J==0 & Q==0) then
                f=cosmo_Schwarzschild;
            elseif (Q==0 & J<>0)
                f=cosmo_Kerr;
            elseif (Q<>0 & J==0)
                f=cosmo_ReissnerNordstrom;
            else
                f=cosmo_KerrNewman;
            end
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                //back to SI units
                Vecc=[Vecc,[Rs/2*ve(3);ve(5);ve(7)]];
                if Conserv==1 then
                    //If Conserv==1 is chosen, then we compute the Hamiltonian and Carter constant at each node
                    g=cosmo_met_mat([0,ve(3),ve(5),ve(7)],rs,rq,a,Lambda);
                    HAM=[HAM,sum([1;ve(4);ve(6);ve(8)].*(g*[1;ve(4);ve(6);ve(8)]))];
                    DT=1+Lambda/3*a^2*cos(ve(5))^2; S=ve(3)^2+a^2*cos(ve(5))^2;
                    Car=S^2*ve(6)^2/DT+chi^2*cos(ve(5))^2/DT*(L0^2/sin(ve(5))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Carter" then
            //Initialize and fix the motion constants
            X=IC; r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6); mu=-Mu; CAR=[];
            Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
            E=-e*r*sqrt(rq)/(chi*S)+sqrt(Dr*Dt*php^2*sin(th)^2/chi^4+(Dr-Dt*a^2*sin(th)^2)/(chi^2*S)*(-mu+rp^2*S/Dr+thp^2*S/Dt));
            Lz=sin(th)^2/(chi^2*(Dt*a^2*sin(th)^2-Dr))*(Dt*(a*chi*(chi*E*(r^2+a^2)+e*r*sqrt(rq))-php*S*Dr)-a*chi^2*E*Dr);
            pt=-E; pr=S*rp/Dr; pth=S*thp/Dt; pph=Lz;
            Q=Dt*pth^2+chi^2*cos(th)^2/Dt*(Lz^2/sin(th)^2-a^2*(E^2+mu*Dt/chi^2+Lambda^2/3*(a*E-Lz)^2));
            k=Q+chi^2*(a*E-Lz)^2;
            X=[r,th,ph,pr,pth]'; Vec=[X];
            //'ode' solve
            function U=F(t,V)
                U=cosmo_Carter_Newman(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                Vecc=[Vecc,[Rs/2*ve(1),ve(2),ve(3)]'];
                if Conserv==1 then
                    gi=cosmo_inv_met_mat([0,ve(1),ve(2),ve(3)],rs,rq,a,Lambda);
                    HAM=[HAM,sum([pt;ve(4);ve(5);pph].*(gi*[pt;ve(4);ve(5);pph]))];
                    DT=1+Lambda/3*a^2*cos(ve(2))^2; S=ve(1)^2+a^2*cos(ve(2))^2;
                    Car=DT*ve(5)^2+chi^2*cos(ve(2))^2/DT*(L0^2/sin(ve(2))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Hamilton" then
            X=cosmo_init_conds_hamiltonian([0,IC(1:3)',time_coeff,IC(4:6)']'); Vec=[X]; CAR=[];
            f=cosmo_Hamilton_equations;
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    gi=cosmo_inv_met_mat([0,ve(2),ve(3),ve(4)],rs,rq,a,Lambda);
                    HAM=[HAM,sum([ve(5);ve(6);ve(7);ve(8)].*(gi*[ve(5);ve(6);ve(7);ve(8)]))];
                    DT=1+Lambda/3*a^2*cos(ve(3))^2; S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=DT*ve(7)^2+chi^2*cos(ve(3))^2/DT*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Symplectic Euler p" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            P=cosmo_met_mat(X,rs,rq,a,Lambda)*P;
            for n=0:N
                //the implicit part is done using 'fsolve'
                function x1=IL(x)
                    x1=x-X-dtau*cosmo_inv_met_mat(x,rs,rq,a,Lambda)*P;
                endfunction
                X=fsolve(X,IL,1e-7); [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(X,rs,rq,a,Lambda);
                P=P-dtau/2*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    DT=1+Lambda/3*a^2*cos(ve(3))^2; S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=DT*Pv(3)^2+chi^2*cos(ve(3))^2/DT*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Symplectic Euler q" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(X,rs,rq,a); P=g*P;
            for n=0:N
                function p1=IR(p)
                    p1=p-P+dtau/2*[0;sum(p.*(drpg*p));sum(p.*(dthpg*p));0]
                endfunction
                P=fsolve(P,IR,1e-7);
                X=X+dtau*cosmo_inv_met_mat(X,rs,rq,a)*P;
                [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(X,rs,rq,a);
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    DT=1+Lambda/3*a^2*cos(ve(3))^2; S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=DT*Pv(3)^2+chi^2*cos(ve(3))^2/DT*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Verlet" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            P=cosmo_met_mat(X,rs,rq,a)*P;
            Ham=[]; [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(X,rs,rq,a);
            for n=[0:N]
                Pp=P-dtau/4*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                X=X+dtau*gi*Pp;
                [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(X,rs,rq,a);
                P=Pp-dtau/4*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    DT=1+Lambda/3*a^2*cos(ve(3))^2; S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=DT*Pv(3)^2+chi^2*cos(ve(3))^2/DT*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Stormer-Verlet" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[0];
            [g,Gam]=cosmo_metric_with_christoffel(X); P=g*P; Ham=[];
            gi=cosmo_inv_met_mat(X,rs,rq,a,Lambda); Ham=[1/2*sum(P.*(gi*P))];
            for n=[1:N]
                function x1=NI(x)
                    x1=x-X-dtau/2*cosmo_inv_met_mat(x,rs,rq,a,Lambda)*P;
                endfunction
                XP=fsolve(X,NI,1e-7); [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(XP,rs,rq,a,Lambda);
                function p1=NJ(p)
                    p1=p-P+dtau/4*([0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0]+[0;sum(p.*(drpg*p));sum(p.*(dthpg*p));0]);
                endfunction
                P=fsolve(P,NJ,1e-7); X=XP+dtau/2*gi*P; Vec=[Vec,X]; Ham=[Ham,sum(P.*(gi*P))]; PVec=[PVec,P];
            end
            Vecc=[]; HAM=[];
            if Conserv==1; then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    DT=1+Lambda/3*a^2*cos(ve(3))^2; S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=DT*Pv(3)^2+chi^2*cos(ve(3))^2/DT*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu*DT/chi^2+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        end
    else
        //The same functions as above, but simplified when Lambda=0 (faster)
        GSI=6.67408e-11; cSI=299792458; e0=8.854187e-12; e=0;
        Mass=cSI^2*Rs/(2*GSI); J=Kerr*GSI*Mass^2/cSI; a=J/(Mass*cSI); 
        Q=Newman*2*Mass*sqrt(%pi*e0*GSI*(1-Kerr^2)); rq2=Q^2*GSI/(4*%pi*e0*cSI^4); rq=4*rq2/Rs^2;
        G=1; c=1; Mass=1; rs=2; rg=1; a=Kerr; Q=sqrt(rq*4*%pi*e0); tau=2*cSI/Rs*Tau;
        IC=[2/Rs*IniConds(1);IniConds(2);IniConds(3);IniConds(4)/cSI;IniConds(5)*Rs/(2*cSI);IniConds(6)*Rs/(2*cSI)];
        dtau=tau/N; r0=IC(1); th0=IC(2); ph0=IC(3); rp0=IC(4); thp0=IC(5); php0=IC(6); S0=r0^2+a^2*cos(th0)^2; D0=r0^2-2*r0+a^2+rq;
        E0=-e*r0*sqrt(rq)/S0+sqrt((D0-a^2*sin(th0)^2)*(S0*rp0^2+S0*D0*thp0^2+D0*Mu)/(S0*D0)+D0*php0^2*sin(th0)^2);
        L0=((S0*D0*php0+a*E0*(rq-2*r0)-a*e*r0*sqrt(rq))*sin(th0)^2)/(D0-a^2*sin(th0)^2);
        time_coeff=((r0^2+a^2)*(E0*(r0^2+a^2)-a*L0+e*r0*sqrt(rq))/D0-a*(a*E0*sin(th0)^2-L0))/S0;
        if Form=="Polar" then
            Vecc=[];
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
                P=atan(y,x)
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
            X=IC; r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6); mu=Mu;
            BB=From_spherical(r,rp,th,thp,ph,php); BBi=[BB(1);BB(3);BB(5)]; BBv=[BB(2);BB(4);BB(6)];
            x=BBi(1); y=BBi(2); z=BBi(3); vx=BBv(1); vy=BBv(2); vz=BBv(3);
            if (abs(z)<1e-10 & abs(vz)<1e-10) then
                th0=0; O0=[1,0,0];
            elseif (abs(z)<1e-10 & abs(vz)>=1e-10) then
                P0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),vz]; Q0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),0];
                th0=sign(vz)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0))); O0=BBi;
            elseif (abs(z)>=1e-10 & abs(vz)<1e-10) then
                th0=%pi/2; O0=BBv;
            else
                O0=[x-z*vx/vz,y-z*vy/vz,0]; P0=[vy*z/vz-y,-vx*z/vz+x,0];
                Q0=[-z*(-vy*z+vz*y)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z*(-vx*z+vz*x)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z];
                th0=sign(z)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0)));
            end
            BBi=rot(O0,-th0,BBi); BBv=rot(O0,-th0,BBv);
            CC=To_spherical(BBi(1),BBv(1),BBi(2),BBv(2),BBi(3),BBv(3));
            R=CC(1); Rp=CC(2); T=CC(3); Tp=CC(4); P=CC(5); Pp=CC(6);
            L=R^2*Pp; X=[1/R;-Rp/L];
            function Y=f(U)
                Y=[U(2);U(1)*(-2*rq*U(1)^2+3*U(1)-1-mu*rq/L^2)+mu/L^2];
            endfunction
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau*Pp:tau*Pp],F);
            Vec=[1./Vec(1,:);P+Pp*[0:dtau:tau]]; Vecc=[]; HAM=[];
            test=0; ii=1;
            while (test==0 & ii<size(Vec)(2))
                Rf=Vec(1,ii); Pf=Vec(2,ii);
                if (Rf<0 | Rf>200*rs) then
                    test=1;
                end
                Cf=rot(O0,th0,[Rf*cos(Pf);Rf*sin(Pf);0]);
                Vecc=[Vecc,[Rs/2*Rf;acos(Cf(3)/Rf);atan(Cf(2),Cf(1))]]; ii=ii+1;
            end
        elseif Form=="Weierstrass" then
            Vecc=[];
            function app=carlson(x,y,z)
                rtol=1e-10; xn=x; yn=y; zn=z; A0=(xn+yn+zn)/3; m=0;
                Q=1/nthroot(3*rtol,6)*max(abs(A0-xn),abs(A0-yn),abs(A0-zn)); A=A0;
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
            X=IC; r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6);
            BB=From_spherical(r,rp,th,thp,ph,php); BBi=[BB(1);BB(3);BB(5)]; BBv=[BB(2);BB(4);BB(6)];
            x=BBi(1); y=BBi(2); z=BBi(3); vx=BBv(1); vy=BBv(2); vz=BBv(3);
            if (abs(z)<1e-10 & abs(vz)<1e-10) then
                th0=0; O0=[1,0,0];
            elseif (abs(z)<1e-10 & abs(vz)>=1e-10) then
                P0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),vz]; Q0=[y*(vx*y-vy*x)/(x^2+y^2),x*(-vx*y+vy*x)/(x^2+y^2),0];
                th0=sign(vz)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0))); O0=BBi;
            elseif (abs(z)>=1e-10 & abs(vz)<1e-10) then
                th0=%pi/2; O0=BBv;
            else
                O0=[x-z*vx/vz,y-z*vy/vz,0]; P0=[vy*z/vz-y,-vx*z/vz+x,0];
                Q0=[-z*(-vy*z+vz*y)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z*(-vx*z+vz*x)*(-vx*y+vy*x)/((vx^2+vy^2)*z^2-2*vz*(vx*x+vy*y)*z+vz^2*(x^2+y^2)),z];
                th0=sign(z)*acos((sum(P0.*Q0))/(norm(P0)*norm(Q0)));
            end
            BBi=rot(O0,-th0,BBi); BBv=rot(O0,-th0,BBv);
            CC=To_spherical(BBi(1),BBv(1),BBi(2),BBv(2),BBi(3),BBv(3));
            R=CC(1); Rp=CC(2); T=CC(3); Tp=CC(4); P=CC(5); Pp=CC(6);
            E=sqrt(Rp^2+(1-2/R+rq/R^2)*mu+Pp^2*(R^2-2*R+rq)); L=R^2*Pp;
            rpol=roots(poly([-rq,2,-1-mu*rq/L^2,2*mu/L^2,(E^2-mu)/L^2],'x','c'));
            frpol=rpol(find(abs(rpol-real(rpol))==min(abs(rpol-real(rpol)))));
            rbar=frpol(find(abs(frpol)==min(abs(frpol)))(1));
            del=(E^2-mu)/L^2; gam=2*(2*del*rbar+mu/L^2); bet=-1-mu*rq/L^2+3*rbar*(gam-2*del*rbar); alp=2+rbar*(2*bet-rbar*(3*gam-4*del*rbar));
            g2=(bet^2/3-alp*gam)/4; g3=(alp*bet*gam/6-alp^2*del/2-bet^3/27)/8;
            function zz=weierP(z)
                N0=12;
                zz0=z/(2^N0); zz=1./zz0^2+g2/20*zz0.^2+g3/28*zz0.^4;
                for j=1:N0
                    zz=-2*zz+(6*zz.^2-g2/2)^2./(4*(4*zz.^3-g2*zz-g3));
                end
            endfunction
            rp2=roots(poly([-g3,-g2,0,4],'u','c')); z0=alp/(4*(R-rbar))+bet/12;
            if abs(Rp)<1e-12 then
                Z0=carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3))
            else
                Z0=sign(-Rp)*carlson(z0-rp2(1),z0-rp2(2),z0-rp2(3))
            end
            deff('r=wrs(t)','r=(2-rbar)*(4*real(weierP(Z0+t))-bet/3)-alp');
            [xxx,vvv,info]=fsolve(0,wrs,1e-15);
            if (abs(vvv)<1e-8 & abs(sign(Pp)*xxx+P)<2*%pi) then
                Z0=-Z0
                tau=min(tau,abs(xxx)); dtau=tau/N;
            end
            Vec=4*real(weierP(Z0+[0:dtau:tau]))-bet/3;
            Vec=[alp./(Vec(1,:))+rbar;P+sign(Pp)*[0:dtau:tau]]; Vecc=[]; HAM=[]; test=0; ii=1;
            while (test==0 & ii<size(Vec)(2))
                Rf=Vec(1,ii); Pf=Vec(2,ii);
                if (Rf<0 | Rf>200*rs) then
                    test=1;
                end
                Cf=rot(O0,th0,[Rf*cos(Pf);Rf*sin(Pf);0]);
                Vecc=[Vecc,[Rs/2*Rf;acos(Cf(3)/Rf);atan(Cf(2),Cf(1))]]; ii=ii+1;
            end
        elseif Form=="Euler-Lagrange" then
            X0=[0,time_coeff,IC(1),IC(4),IC(2),IC(5),IC(3),IC(6)]'; Vec=[X0]; X=X0; CAR=[];
            if (J==0 & Q==0) then
                f=Schwarzschild;
            elseif (Q==0 & J<>0)
                f=Kerr;
            elseif (Q<>0 & J==0)
                f=ReissnerNordstrom;
            else
                f=KerrNewman;
            end
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                Vecc=[Vecc,[Rs/2*ve(3);ve(5);ve(7)]];
                if Conserv==1 then
                    g=met_mat([0,ve(3),ve(5),ve(7)],rs,rq,a);
                    HAM=[HAM,sum([1;ve(4);ve(6);ve(8)].*(g*[1;ve(4);ve(6);ve(8)]))];
                    S=ve(3)^2+a^2*cos(ve(5))^2;
                    Car=S^2*ve(6)^2+cos(ve(5))^2*(L0^2/sin(ve(5))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Carter" then
            e=0;
            X=IC; r=X(1); th=X(2); ph=X(3); rp=X(4); thp=X(5); php=X(6); mu=-Mu; CAR=[];
            S=r^2+a^2*cos(th)^2; D=r^2-2*r+a^2+rq;
            E=-e*r*sqrt(rq)/S+sqrt(D*php^2*sin(th)^2+(a^2*sin(th)^2-D)*(mu-rp^2*S/D-S*thp^2)/S);
            Lz=((S*D*php+a*E*(rq-2*r)-a*e*r*sqrt(rq))*sin(th)^2)/(D-a^2*sin(th)^2);
            pt=-E; pr=S*rp/D; pth=S*thp; pph=Lz;
            Q=pth^2+cos(th)^2*(Lz^2/sin(th)^2-a^2*(E^2+mu));
            k=Q+Lz^2+a^2*(E^2+mu);
            X=[r,th,ph,pr,pth]'; Vec=[X];
            f=Carter_Newman;
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                Vecc=[Vecc,[Rs/2*ve(1),ve(2),ve(3)]'];
                if Conserv==1 then
                    gi=inv_met_mat([0,ve(1),ve(2),ve(3)],rs,rq,a);
                    HAM=[HAM,sum([pt;ve(4);ve(5);pph].*(gi*[pt;ve(4);ve(5);pph]))];
                    S=ve(1)^2+a^2*cos(ve(2))^2;
                    Car=ve(5)^2+cos(ve(2))^2*(L0^2/sin(ve(2))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Hamilton" then
            X=init_conds_hamiltonian([0,IC(1:3)',time_coeff,IC(4:6)']'); Vec=[X]; CAR=[];
            f=Hamilton_equations;
            function U=F(t,V)
                U=f(V);
            endfunction
            Vec=ode("adams",X,0,[0:dtau:tau],F);
            Vecc=[]; HAM=[];
            for ve=Vec
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    gi=inv_met_mat([0,ve(2),ve(3),ve(4)],rs,rq,a);
                    HAM=[HAM,sum([ve(5);ve(6);ve(7);ve(8)].*(gi*[ve(5);ve(6);ve(7);ve(8)]))];
                    S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=ve(7)^2+cos(ve(3))^2*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Symplectic Euler p" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[P]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            P=met_mat(X,rs,rq,a)*P;
            for n=0:N
                function x1=IL(x)
                    x1=x-X-dtau*inv_met_mat(x,rs,rq,a)*P;
                endfunction
                X=fsolve(X,IL,1e-7); [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a);
                P=P-dtau/2*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=Pv(3)^2+cos(ve(3))^2*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Symplectic Euler q" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[P]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a); P=g*P;
            for n=0:N
                function p1=IR(p)
                    p1=p-P+dtau/2*[0;sum(p.*(drpg*p));sum(p.*(dthpg*p));0]
                endfunction
                P=fsolve(P,IR,1e-7);
                X=X+dtau*inv_met_mat(X,rs,rq,a)*P;
                [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a);
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=Pv(3)^2+cos(ve(3))^2*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Verlet" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[P]; Ham=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[];
            P=met_mat(X,rs,rq,a)*P;
            Ham=[]; [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a);
            for n=[0:N]
                Pp=P-dtau/4*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                X=X+dtau*gi*Pp;
                [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a);
                P=Pp-dtau/4*[0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0];
                Vec=[Vec,X]; PVec=[PVec,P]; Ham=[Ham;sum(P.*(g*P))];
            end
            Vecc=[]; HAM=[];
            if Conserv==1 then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=Pv(3)^2+cos(ve(3))^2*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        elseif Form=="Stormer-Verlet" then
            P=[time_coeff,IC(4:6)']'; X=[0,IC(1:3)']'; Vec=[]; PVec=[]; k=Kerr; rq=4*rq2/Rs^2; rs=2; CAR=[0];
            [g,gi,drpg,dthpg]=inverse_metric_matrix(X,rs,rq,a); P=g*P; 
            Ham=[1/2*sum(P.*(gi*P))];
            for n=[1:N]
                function x1=NI(x)
                    x1=x-X-dtau/2*inv_met_mat(x,rs,rq,a)*P;
                endfunction
                XP=fsolve(X,NI,1e-7); [g,gi,drpg,dthpg]=inverse_metric_matrix(XP,rs,rq,a);
                function p1=NJ(p)
                    p1=p-P+dtau/4*([0;sum(P.*(drpg*P));sum(P.*(dthpg*P));0]+[0;sum(p.*(drpg*p));sum(p.*(dthpg*p));0]);
                endfunction
                P=fsolve(P,NJ,1e-7); X=XP+dtau/2*gi*P; Vec=[Vec,X]; Ham=[Ham,sum(P.*(gi*P))]; PVec=[PVec,P];
            end
            Vecc=[]; HAM=[];
            if Conserv==1; then
                HAM=Ham;
            end
            for vee=1:size(Vec)(2)
                ve=Vec(:,vee);
                Vecc=[Vecc,[Rs/2*ve(2);ve(3);ve(4)]];
                if Conserv==1 then
                    Pv=PVec(:,vee);
                    S=ve(2)^2+a^2*cos(ve(3))^2;
                    Car=Pv(3)^2+cos(ve(3))^2*(L0^2/sin(ve(3))^2-a^2*(E0^2-Mu+Lambda^2/3*(a*E0-L0)^2));
                    CAR=[CAR,Car];
                end
            end
        end
    end
endfunction

