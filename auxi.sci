//Fundamental constants
G=1; c=1; cSI=299792458; GSI=6.67408e-11; e0=8.854187e-12;

//Compute (x,y,z) and its derivative (xp,yp,zp) from Boyer-Lidquist coordinates (R,T,P) and (Rp,Tp,Pp).
//We make the dependance in the Kerr parameter explicit.
function BB=BoyerLindquist(R,Rp,T,Tp,P,Pp,a)
    x=sqrt(R^2+a^2)*sin(T)*cos(P);
    y=sqrt(R^2+a^2)*sin(T)*sin(P);
    z=R*cos(T);
    xp=(Tp*cos(T)*cos(P)*(R^2+a^2)-sin(T)*(Pp*sin(P)*(R^2+a^2)-R*Rp*cos(P)))/sqrt(R^2+a^2);
    yp=(Tp*cos(T)*sin(P)*(R^2+a^2)+sin(T)*(Pp*cos(P)*(R^2+a^2)+R*Rp*sin(P)))/sqrt(R^2+a^2);
    zp=Rp*cos(T)-R*Tp*sin(T);
    BB=[x,xp,y,yp,z,zp]';
endfunction

//The inverse transformation, giving (R,T,P) and (Rp,Tp,Pp) as a function of (x,y,z) and (xp,yp,zp)
function XX=InvBoyerLindquist(x,xp,y,yp,z,zp,a)
    P=atan(y,x);
    R=sqrt((-a^2+x^2+y^2+z^2+sqrt(a^2*(a^2-2*x^2-2*y^2+2*z^2)+(x^2+y^2+z^2)^2))/2);
    T=acos(z/R);
    Rp=R*(x*xp+y*yp+z*zp)/(2*R^2+a^2-x^2-y^2-z^2)+a^2*z*zp/(R*(2*R^2+a^2-x^2-y^2-z^2));
    Tp=(z*Rp-zp*R)/(R*sqrt(R^2-z^2));
    Pp=(yp*x-xp*y)/(x^2+y^2);
    XX=[R,Rp,T,Tp,P,Pp]';
endfunction







//The metric matrix depending on the three parameters and a vector V=(t,r,theta,phi)
function g=met_mat(V,rs,rq,a)
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & a<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r-rq)/(r^2+a^2*cos(theta)^2), 0, 0, (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r+rq), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(r^2+2*r-rq)*a^2+r^4)/(r^2+a^2*cos(theta)^2)];
    elseif (rq==0 & a<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r)/(r^2+a^2*cos(theta)^2), 0, 0, -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r)*cos(theta)^2+((r+2)*a^2+r^3)*r)/(r^2+a^2*cos(theta)^2)];
    elseif (rq<>0 & a==0) then
        g=diag([(-r^2+2*r-rq)/r^2, r^2/(r^2-2*r+rq), r^2, sin(theta)^2*r^2]);
    else
        g=diag([(-r+2)/r, r/(r-2), r^2, sin(theta)^2*r^2]);
    end
endfunction

//The metric matrix, its r and theta derivatives and its inverse depending on the three parameters and a vector V=(t,r,theta,phi)
function [g,gi,drg,dthg]=metric_matrix(V,rs,rq,a)
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & a<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r-rq)/(r^2+a^2*cos(theta)^2), 0, 0, (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r+rq), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(r^2+2*r-rq)*a^2+r^4)/(r^2+a^2*cos(theta)^2)];
        gi=[(-a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(-r^2-2*r+rq)*a^2-r^4)/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r+rq)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r+rq)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r+rq)*sin(theta)^2)];
        drg=[(2*a^2*cos(theta)^2-2*r^2+2*r*rq)/(r^2+a^2*cos(theta)^2)^2, 0, 0, -2*a*sin(theta)^2*(a^2*cos(theta)^2-r^2+r*rq)/(r^2+a^2*cos(theta)^2)^2; 0, (-2*a^2*(r-1)*cos(theta)^2+2*r*(a^2-r+rq))/(a^2+r^2-2*r+rq)^2, 0, 0; 0, 0, 2*r, 0; -2*a*sin(theta)^2*(a^2*cos(theta)^2-r^2+r*rq)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 2*sin(theta)^2*(a^4*(r-1)*cos(theta)^4+a^2*(2*r^3+a^2+r^2-r*rq)*cos(theta)^2-((r-rq)*a^2-r^4)*r)/(r^2+a^2*cos(theta)^2)^2];
        dthg=[(4*r-2*rq)*a^2*sin(theta)*cos(theta)/(r^2+a^2*cos(theta)^2)^2, 0, 0, -4*cos(theta)*(a^2+r^2)*a*sin(theta)*(r-(1/2)*rq)/(r^2+a^2*cos(theta)^2)^2; 0, -2*a^2*cos(theta)*sin(theta)/(a^2+r^2-2*r+rq), 0, 0; 0, 0, -2*a^2*cos(theta)*sin(theta), 0; -4*cos(theta)*(a^2+r^2)*a*sin(theta)*(r-(1/2)*rq)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 2*cos(theta)*sin(theta)*(a^4*(a^2+r^2-2*r+rq)*cos(theta)^4+2*a^2*r^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(2*r-rq)*a^4+r^2*(r^2+4*r-2*rq)*a^2+r^6)/(r^2+a^2*cos(theta)^2)^2];
    elseif (rq==0 & a<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r)/(r^2+a^2*cos(theta)^2), 0, 0, -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r)*cos(theta)^2+((r+2)*a^2+r^3)*r)/(r^2+a^2*cos(theta)^2)];
        gi=[(-a^2*(a^2+r^2-2*r)*cos(theta)^2-a^2*r^2-r^4-2*a^2*r)/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r)*sin(theta)^2)];
        drg=[(2*a^2*cos(theta)^2-2*r^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0, -2*a*sin(theta)^2*(a^2*cos(theta)^2-r^2)/(r^2+a^2*cos(theta)^2)^2; 0, (-2*a^2*(r-1)*cos(theta)^2+2*a^2*r-2*r^2)/(a^2+r^2-2*r)^2, 0, 0; 0, 0, 2*r, 0; -2*a*sin(theta)^2*(a^2*cos(theta)^2-r^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 2*sin(theta)^2*(a^4*(r-1)*cos(theta)^4+a^2*(2*r^3+a^2+r^2)*cos(theta)^2+r^5-a^2*r^2)/(r^2+a^2*cos(theta)^2)^2];
        dthg=[4*a^2*cos(theta)*sin(theta)*r/(r^2+a^2*cos(theta)^2)^2, 0, 0, -4*r*a*sin(theta)*cos(theta)*(a^2+r^2)/(r^2+a^2*cos(theta)^2)^2; 0, -2*a^2*cos(theta)*sin(theta)/(a^2+r^2-2*r), 0, 0; 0, 0, -2*a^2*cos(theta)*sin(theta), 0; -4*r*a*sin(theta)*cos(theta)*(a^2+r^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 2*sin(theta)*(a^4*(a^2+r^2-2*r)*cos(theta)^4+2*a^2*r^2*(a^2+r^2-2*r)*cos(theta)^2+2*a^4*r+r^3*(r+4)*a^2+r^6)*cos(theta)/(r^2+a^2*cos(theta)^2)^2];
    elseif (rq<>0 & a==0) then
        g=[(-r^2+2*r-rq)/r^2, 0, 0, 0; 0, r^2/(r^2-2*r+rq), 0, 0; 0, 0, r^2, 0; 0, 0, 0, sin(theta)^2*r^2];
        gi=[-r^2/(r^2-2*r+rq), 0, 0, 0; 0, (r^2-2*r+rq)/r^2, 0, 0; 0, 0, 1/r^2, 0; 0, 0, 0, 1/(r^2*sin(theta)^2)];
        drg=[(-2*r+2*rq)/r^3, 0, 0, 0; 0, -2*r*(r-rq)/(r^2-2*r+rq)^2, 0, 0; 0, 0, 2*r, 0; 0, 0, 0, 2*sin(theta)^2*r];
        dthg=[0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 2*sin(theta)*r^2*cos(theta)];
    else
        g=[(-r+2)/r, 0, 0, 0; 0, r/(r-2), 0, 0; 0, 0, r^2, 0; 0, 0, 0, sin(theta)^2*r^2];
        gi=[-r/(r-2), 0, 0, 0; 0, (r-2)/r, 0, 0; 0, 0, 1/r^2, 0; 0, 0, 0, 1/(r^2*sin(theta)^2)];
        drg=[-2/r^2, 0, 0, 0; 0, -2/(r-2)^2, 0, 0; 0, 0, 2*r, 0; 0, 0, 0, 2*sin(theta)^2*r];
        dthg=[0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 2*sin(theta)*r^2*cos(theta)];
    end
endfunction

//Inverse metric matrix
function gi=inv_met_mat(V,rs,rq,a)
    r=V(2); theta=V(3); phi=V(4);
    if (a<>0 & rq<>0) then
        gi=[(-a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(-r^2-2*r+rq)*a^2-r^4)/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r+rq)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r+rq)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r+rq)*sin(theta)^2)];
    elseif (a<>0 & rq==0)
        gi=[(-a^2*(a^2+r^2-2*r)*cos(theta)^2-a^2*r^2-r^4-2*a^2*r)/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r)*sin(theta)^2)];
    elseif (a==0 & rq<>0)
        gi=[-r^2/(r^2-2*r+rq), 0, 0, 0; 0, (r^2-2*r+rq)/r^2, 0, 0; 0, 0, 1/r^2, 0; 0, 0, 0, 1/(r^2*sin(theta)^2)];
    else
        gi=[-r/(r-2), 0, 0, 0; 0, (r-2)/r, 0, 0; 0, 0, 1/r^2, 0; 0, 0, 0, 1/(r^2*sin(theta)^2)];
    end
endfunction

//Inverse metric matrix and its r and theta derivatives as a function of the three parameters
function [g,gi,drpg,dthpg]=inverse_metric_matrix(V,rs,rq,a)
    r=V(2); theta=V(3); phi=V(4);
    if (a<>0 & rq<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r-rq)/(r^2+a^2*cos(theta)^2), 0, 0, (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r+rq), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(r^2+2*r-rq)*a^2+r^4)/(r^2+a^2*cos(theta)^2)];
        gi=[(-a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(-r^2-2*r+rq)*a^2-r^4)/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r+rq)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -(2*r-rq)*a/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r+rq)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r+rq)*sin(theta)^2)];
        drpg=[(-2*(r^4-4*r^3+(2*a^2+4*rq)*r^2-r*rq^2+a^4)*a^2*cos(theta)^2+2*(r^5-r^4*rq+2*a^2*r^3-2*a^2*(rq+2)*r^2+a^2*(a^2+4*rq)*r-a^4*rq-a^2*rq^2)*r)/((a^2+r^2-2*r+rq)^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -2*a*(a^2*(a^2-r^2+r*rq)*cos(theta)^2-(3*r^3+(-2*rq-4)*r^2+(a^2+4*rq)*r-a^2*rq-rq^2)*r)/((a^2+r^2-2*r+rq)^2*(r^2+a^2*cos(theta)^2)^2); 0, (2*a^2*(r-1)*cos(theta)^2-2*r*(a^2-r+rq))/(r^2+a^2*cos(theta)^2)^2, 0, 0; 0, 0, -2*r/(r^2+a^2*cos(theta)^2)^2, 0; -2*a*(a^2*(a^2-r^2+r*rq)*cos(theta)^2-(3*r^3+(-2*rq-4)*r^2+(a^2+4*rq)*r-a^2*rq-rq^2)*r)/((a^2+r^2-2*r+rq)^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2*a^4*(r-1)*cos(theta)^4+2*a^2*(2*r^3+a^2-3*r^2+r*rq)*cos(theta)^2-2*(-r^4+4*r^3+(-2*rq-4)*r^2+(a^2+4*rq)*r-a^2*rq-rq^2)*r)/((a^2+r^2-2*r+rq)^2*(r^2+a^2*cos(theta)^2)^2*(cos(theta)^2-1))];
        dthpg=[-4*cos(theta)*(a^2+r^2)*a^2*(r-(1/2)*rq)*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq)), 0, 0, -2*a^3*cos(theta)*sin(theta)*(2*r-rq)/((r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq)); 0, 2*(a^2+r^2-2*r+rq)*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2)^2, 0, 0; 0, 0, 2*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2)^2, 0; -2*a^3*cos(theta)*sin(theta)*(2*r-rq)/((r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq)), 0, 0, 2*cos(theta)*(cos(theta)^4*a^4+2*a^2*(r^2-2*r+rq)*cos(theta)^2+r^4+2*a^2*r-a^2*rq-2*r^3+r^2*rq)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq)*(cos(theta)^2-1))];
    elseif (a<>0 & rq==0)
        g=[(-a^2*cos(theta)^2-r^2+2*r)/(r^2+a^2*cos(theta)^2), 0, 0, -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r)*cos(theta)^2+((r+2)*a^2+r^3)*r)/(r^2+a^2*cos(theta)^2)];
        gi=[(-a^2*(a^2+r^2-2*r)*cos(theta)^2-a^2*r^2-r^4-2*a^2*r)/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)); 0, (a^2+r^2-2*r)/(r^2+a^2*cos(theta)^2), 0, 0; 0, 0, 1/(r^2+a^2*cos(theta)^2), 0; -2*a*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)), 0, 0, (a^2*cos(theta)^2+r^2-2*r)/((r^2+a^2*cos(theta)^2)*(a^2+r^2-2*r)*sin(theta)^2)];
        drpg=[(-2*a^2*(a^4+2*a^2*r^2+r^4-4*r^3)*cos(theta)^2+2*a^4*r^2+4*a^2*r^4+2*r^6-8*a^2*r^3)/((a^2+r^2-2*r)^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -2*a*((a^4-a^2*r^2)*cos(theta)^2-a^2*r^2-3*r^4+4*r^3)/((a^2+r^2-2*r)^2*(r^2+a^2*cos(theta)^2)^2); 0, (2*a^2*(r-1)*cos(theta)^2-2*a^2*r+2*r^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0; 0, 0, -2*r/(r^2+a^2*cos(theta)^2)^2, 0; -2*a*((a^4-a^2*r^2)*cos(theta)^2-a^2*r^2-3*r^4+4*r^3)/((a^2+r^2-2*r)^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2*a^4*(r-1)*cos(theta)^4+2*a^2*(2*r^3+a^2-3*r^2)*cos(theta)^2-2*r^2*(-r^3+a^2+4*r^2-4*r))/((a^2+r^2-2*r)^2*(r^2+a^2*cos(theta)^2)^2*(cos(theta)^2-1))];
        dthpg=[-4*sin(theta)*r*cos(theta)*a^2*(a^2+r^2)/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2), 0, 0, -4*cos(theta)*sin(theta)*a^3*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2); 0, 2*(a^2+r^2-2*r)*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2)^2, 0, 0; 0, 0, 2*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2)^2, 0; -4*cos(theta)*sin(theta)*a^3*r/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2), 0, 0, 2*(cos(theta)^4*a^4+2*a^2*r*(r-2)*cos(theta)^2+r^4+2*a^2*r-2*r^3)*cos(theta)/(sin(theta)*(a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2*(cos(theta)^2-1))];
    elseif (a==0 & rq<>0)
        g=diag([-c^2*(r^2-r*rs+rq)/r^2,r^2/(r^2-r*rs+rq),r^2,r^2*sin(theta)^2]);
        gi=diag([-r^2/(c^2*(r^2-r*rs+rq)),(r^2-r*rs+rq)/r^2,1/r^2,1/(r^2*sin(theta)^2)]);
        drpg=diag([r*(r*rs-2*rq)/(c^2*(r^2-r*rs+rq)^2),(r*rs-2*rq)/r^3,-2/r^3,-2/(r^3*sin(theta)^2)]);
        dthpg=diag([0, 0, 0, -2*cos(theta)/(r^2*sin(theta)^3)]);
    else
        g=diag([-c^2*(1-rs/r),1/(1-rs/r),r^2,r^2*sin(theta)^2]);
        gi=diag([-1/(c^2*(1-rs/r)),1-rs/r,1/r^2,1/(r^2*sin(theta)^2)]);
        drpg=diag([rs/(c^2*(r-rs)^2),rs/r^2,-2/r^3,-2/(r^3*sin(theta)^2)]);
        dthpg=diag([0,0,0,-2*cos(theta)/(r^2*sin(theta)^3)]);
    end
endfunction

//Metric matrix with Christoffel symbols
function [g,Gam]=metric_with_christoffel(V)
    G=1; c=1;
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & k<>0 & rs<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r-rq)/(r^2+a^2*cos(theta)^2), 0, 0, (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r+rq), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; (-2*r+rq)*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(r^2+2*r-rq)*a^2+r^4)/(r^2+a^2*cos(theta)^2)];
        Gam=list([0, -2*(a^2*cos(theta)^2-r^2+r*rq)*(a^2+r^2)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r+rq)), (-4*r+2*rq)*a^2*sin(theta)*cos(theta)/(a^2*cos(theta)^2+r^2)^2, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 2*sin(theta)^2*a*(a^2*(a^2-r^2+r*rq)*cos(theta)^2-a^2*r^2+a^2*r*rq-3*r^4+2*r^3*rq)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r+rq)), (4*r-2*rq)*a^3*sin(theta)^3*cos(theta)/(a^2*cos(theta)^2+r^2)^2, 0], [-(a^2+r^2-2*r+rq)*(a^2*cos(theta)^2-r^2+r*rq)/(a^2*cos(theta)^2+r^2)^3, 0, 0, (a^2+r^2-2*r+rq)*a*sin(theta)^2*(a^2*cos(theta)^2-r^2+r*rq)/(a^2*cos(theta)^2+r^2)^3; 0, (-a^2*(r-1)*cos(theta)^2+r*(a^2-r+rq))/((a^2+r^2-2*r+rq)*(a^2*cos(theta)^2+r^2)), -2*a^2*cos(theta)*sin(theta)/(a^2*cos(theta)^2+r^2), 0; 0, 0, -(a^2+r^2-2*r+rq)*r/(a^2*cos(theta)^2+r^2), 0; (a^2+r^2-2*r+rq)*a*sin(theta)^2*(a^2*cos(theta)^2-r^2+r*rq)/(a^2*cos(theta)^2+r^2)^3, 0, 0, -sin(theta)^2*(a^2+r^2-2*r+rq)*(a^4*(r-1)*cos(theta)^4+a^2*(2*r^3+a^2+r^2-r*rq)*cos(theta)^2-((r-rq)*a^2-r^4)*r)/(a^2*cos(theta)^2+r^2)^3], [(-2*r+rq)*a^2*sin(theta)*cos(theta)/(a^2*cos(theta)^2+r^2)^3, 0, 0, (2*r-rq)*a*sin(theta)*cos(theta)*(a^2+r^2)/(a^2*cos(theta)^2+r^2)^3; 0, a^2*cos(theta)*sin(theta)/((a^2+r^2-2*r+rq)*(a^2*cos(theta)^2+r^2)), 0, 0; 0, 2*r/(a^2*cos(theta)^2+r^2), -a^2*cos(theta)*sin(theta)/(a^2*cos(theta)^2+r^2), 0; (2*r-rq)*a*sin(theta)*cos(theta)*(a^2+r^2)/(a^2*cos(theta)^2+r^2)^3, 0, 0, -cos(theta)*sin(theta)*(a^4*(a^2+r^2-2*r+rq)*cos(theta)^4+2*a^2*r^2*(a^2+r^2-2*r+rq)*cos(theta)^2+(2*r-rq)*a^4+r^2*(r^2+4*r-2*rq)*a^2+r^6)/(a^2*cos(theta)^2+r^2)^3], [0, -2*a*(a^2*cos(theta)^2-r^2+r*rq)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r+rq)), (-4*r+2*rq)*a*cos(theta)/((a^2*cos(theta)^2+r^2)^2*sin(theta)), 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, (2*a^4*(r-1)*cos(theta)^4+2*a^2*(2*r^3+a^2-r^2)*cos(theta)^2+2*r^5-2*a^2*r^2+2*a^2*r*rq-4*r^4+2*r^3*rq)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r+rq)), 2*(cos(theta)^4*a^4+2*a^2*(r^2-r+(1/2)*rq)*cos(theta)^2+a^2*(2*r-rq)+r^4)*cos(theta)/((a^2*cos(theta)^2+r^2)^2*sin(theta)), 0]);
    elseif (rq==0 & k<>0 & rs<>0) then
        g=[(-a^2*cos(theta)^2-r^2+2*r)/(r^2+a^2*cos(theta)^2), 0, 0, -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2); 0, (r^2+a^2*cos(theta)^2)/(a^2+r^2-2*r), 0, 0; 0, 0, r^2+a^2*cos(theta)^2, 0; -2*r*a*sin(theta)^2/(r^2+a^2*cos(theta)^2), 0, 0, sin(theta)^2*(a^2*(a^2+r^2-2*r)*cos(theta)^2+((r+2)*a^2+r^3)*r)/(r^2+a^2*cos(theta)^2)];
        Gam=list([0, -2*(a^2+r^2)*(a^2*cos(theta)^2-r^2)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r)), -4*a^2*cos(theta)*sin(theta)*r/(a^2*cos(theta)^2+r^2)^2, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 2*((a^4-a^2*r^2)*cos(theta)^2-a^2*r^2-3*r^4)*a*sin(theta)^2/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r)), 4*r*a^3*sin(theta)^3*cos(theta)/(a^2*cos(theta)^2+r^2)^2, 0], [-(a^2+r^2-2*r)*(a^2*cos(theta)^2-r^2)/(a^2*cos(theta)^2+r^2)^3, 0, 0, (a^2+r^2-2*r)*a*sin(theta)^2*(a^2*cos(theta)^2-r^2)/(a^2*cos(theta)^2+r^2)^3; 0, (-a^2*(r-1)*cos(theta)^2+r*(a^2-r))/((a^2+r^2-2*r)*(a^2*cos(theta)^2+r^2)), -2*a^2*cos(theta)*sin(theta)/(a^2*cos(theta)^2+r^2), 0; 0, 0, -(a^2+r^2-2*r)*r/(a^2*cos(theta)^2+r^2), 0; (a^2+r^2-2*r)*a*sin(theta)^2*(a^2*cos(theta)^2-r^2)/(a^2*cos(theta)^2+r^2)^3, 0, 0, -sin(theta)^2*(a^2+r^2-2*r)*(a^4*(r-1)*cos(theta)^4+a^2*(2*r^3+a^2+r^2)*cos(theta)^2+r^5-a^2*r^2)/(a^2*cos(theta)^2+r^2)^3], [-2*a^2*cos(theta)*sin(theta)*r/(a^2*cos(theta)^2+r^2)^3, 0, 0, 2*r*a*sin(theta)*cos(theta)*(a^2+r^2)/(a^2*cos(theta)^2+r^2)^3; 0, a^2*cos(theta)*sin(theta)/((a^2*cos(theta)^2+r^2)*(a^2+r^2-2*r)), 0, 0; 0, 2*r/(a^2*cos(theta)^2+r^2), -a^2*cos(theta)*sin(theta)/(a^2*cos(theta)^2+r^2), 0; 2*r*a*sin(theta)*cos(theta)*(a^2+r^2)/(a^2*cos(theta)^2+r^2)^3, 0, 0, -(a^4*(a^2+r^2-2*r)*cos(theta)^4+2*a^2*r^2*(a^2+r^2-2*r)*cos(theta)^2+2*a^4*r+r^3*(r+4)*a^2+r^6)*sin(theta)*cos(theta)/(a^2*cos(theta)^2+r^2)^3], [0, (-2*cos(theta)^2*a^3+2*a*r^2)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r)), -4*cos(theta)*a*r/(sin(theta)*(a^2*cos(theta)^2+r^2)^2), 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, (2*a^4*(r-1)*cos(theta)^4+2*a^2*(2*r^3+a^2-r^2)*cos(theta)^2+2*r^5-2*a^2*r^2-4*r^4)/((a^2*cos(theta)^2+r^2)^2*(a^2+r^2-2*r)), 2*(cos(theta)^4*a^4+2*a^2*r*(r-1)*cos(theta)^2+r^4+2*a^2*r)*cos(theta)/(sin(theta)*(a^2*cos(theta)^2+r^2)^2), 0]);
    elseif (rq<>0 & k==0 & rs<>0) then
        g=diag([-(r^2-2*r+rq)/r^2, r^2/(r^2-2*r+rq), r^2, r^2*sin(theta)^2]);
        Gam=list([0, (2*r - 2*rq)/(r*(r^2 - 2*r + rq)), 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0], [(r^2 - 2*r + rq)*(r - rq)/r^5, 0, 0, 0; 0, (-r + rq)/(r*(r^2 - 2*r + rq)), 0, 0; 0, 0, (-r^2 + 2*r - rq)/r, 0; 0, 0, 0, -(r^2 - 2*r + rq)*sin(theta)^2/r], [0, 0, 0, 0; 0, 0, 0, 0; 0, 2/r, 0, 0; 0, 0, 0, -cos(theta)*sin(theta)], [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 2/r, 2*cos(theta)/sin(theta), 0]);
    elseif (rq==0 & k==0 & rs<>0) then
        g=diag([-(r^2-2*r)/r^2, r^2/(r^2-2*r), r^2, r^2*sin(theta)^2]);
        Gam=list([0, 2/(r*(r-2)), 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0], [(r-2)/r^3, 0, 0, 0; 0, -1/(r*(r-2)), 0, 0; 0, 0, -r+2, 0; 0, 0, 0, -(r-2)*sin(theta)^2], [0, 0, 0, 0; 0, 0, 0, 0; 0, 2/r, 0, 0; 0, 0, 0, -cos(theta)*sin(theta)], [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 2/r, 2*cos(theta)/sin(theta), 0]);
    elseif rs==0 then
        g=diag([-1,1,1,1]); Gam=list(zeros(4,4),zeros(4,4),zeros(4,4),zeros(4,4));
    end
endfunction











//Hamilton equations
//The initial conditions for Hamilton equations
function X=init_conds_hamiltonian(V)
    g=met_mat(V,rs,rq,a);
    X=[V(1:4)',(g*V(5:8))']';
endfunction

//The function defining the differential system of Hamilton's equations
function Y=Hamilton_equations(V)
    [g,gi,drpg,dthpg]=inverse_metric_matrix(V(1:4),rs,rq,a);
    Y=[gi(1,:)*V(5:8),gi(2,:)*V(5:8),gi(3,:)*V(5:8),gi(4,:)*V(5:8),0,-sum(V(5:8).*(drpg*V(5:8)))/2,-sum(V(5:8).*(dthpg*V(5:8)))/2,0]';
endfunction











//Carter-Newman (inspired by PYYY)
//The function defining the Carter's equations
function Y=Carter_Newman(V)
    ep=1e-10;
    r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
    D=r^2-2*r+a^2+rq; S=r^2+a^2*cos(th)^2;
    if abs(D)<ep then
        D=sign(D)*ep;
    end
    if abs(S)<ep then
        S=sign(S)*ep;
    end
    if abs(sin(th))<ep then
        th=sign(th)*asin(ep);
    end
    prp=((r-1)*(mu*(r^2+a^2)-k)+r*D*mu+2*r*E^2*(r^2+a^2)-2*a*E*Lz)/(S*D)-(2*pr^2*(r-1))/S + (e*r*sqrt(rq)+(a^2+3*r^2)*E-a*Lz)*e*sqrt(rq)/(S*D);
    pthp=sin(th)*cos(th)*(Lz^2/sin(th)^4-a^2*(E^2+mu))/S;
    rp=D*pr/S;
    thp=pth/S;
    php=a*(E*(r^2+a^2)-a*Lz+e*r*sqrt(rq))/(S*D)+(Lz/sin(th)^2-a*E)/S;
    Y=[rp,thp,php,prp,pthp]';
endfunction














//Euler-Lagrange equation (computed with simple Christoffel symbols in Maple)

//The function defining the differential system of the Euler-Lagrange equations
//For Kerr-Newman (a<>0<>Q)
function Y=KerrNewman(V)//Q<>0, J<>0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, (2*a^3*r1*phi1*(a^2-r*(r-rq))*cos(theta)^4+4*sin(theta)*a^3*(r-(1/2)*rq)*(a^2+r^2-2*r+rq)*phi1*theta1*cos(theta)^3-2*(a^4*phi1-a^3*t1-a*r^2*t1+3*(r-(2/3)*rq)*phi1*r^3)*a*r1*cos(theta)^2-4*sin(theta)*(a*phi1-t1)*a^2*(r-(1/2)*rq)*(a^2+r^2-2*r+rq)*theta1*cos(theta)+2*r1*(phi1*(r-rq)*a^3-t1*(r-rq)*a^2+3*(r-(2/3)*rq)*phi1*r^2*a-r^2*t1*(r-rq))*r)/((r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq)), r1, (-a^4*(r-1)*(a^2*phi1+phi1*r^2-a*r1-2*phi1*r+phi1*rq)*(a^2*phi1+phi1*r^2+a*r1-2*phi1*r+phi1*rq)*cos(theta)^6+2*a^6*sin(theta)*r1*theta1*(a^2+r^2-2*r+rq)*cos(theta)^5+(-2*phi1^2*r^7+7*phi1^2*r^6+((-3*phi1^2+theta1^2)*a^2-3*(rq+4/3)*phi1^2)*r^5+(-4*a^2*theta1^2+2*a*t1*phi1+2*phi1^2*(rq-2))*r^4+(2*a^4*theta1^2+(2*rq*theta1^2+16*phi1^2+2*r1^2+4*theta1^2)*a^2-8*a*t1*phi1+8*rq*phi1^2)*r^3+((-9*phi1^2-4*theta1^2)*a^4+4*a^3*t1*phi1+((-14*rq-8)*phi1^2-4*rq*theta1^2-r1^2)*a^2+4*t1*phi1*(rq+2)*a-5*rq^2*phi1^2)*r^2+(a^2+rq)*((phi1^2+theta1^2)*a^4+((2*rq+8)*phi1^2+rq*theta1^2-r1^2)*a^2-8*a*t1*phi1+rq^2*phi1^2)*r-2*a*phi1*(a^2+rq)^2*(a*phi1-t1))*a^2*cos(theta)^4+4*a^4*sin(theta)*r^2*r1*theta1*(a^2+r^2-2*r+rq)*cos(theta)^3+(-phi1^2*r^9+4*phi1^2*r^8+(2*a^2*theta1^2-2*phi1^2*(rq+2))*r^7+((-2*phi1^2-8*theta1^2)*a^2-2*a*t1*phi1+4*rq*phi1^2)*r^6+((3*phi1^2+4*theta1^2)*a^4+(4*rq*theta1^2+r1^2+8*theta1^2)*a^2+2*t1*phi1*(rq+4)*a-rq^2*phi1^2)*r^5-3*a*((phi1^2+(8/3)*theta1^2)*a^3+2*a^2*t1*phi1+((-(4/3)*rq-8/3)*phi1^2+(8/3)*rq*theta1^2-(1/3)*r1^2-(1/3)*t1^2)*a+4*t1*(rq+2/3)*phi1)*r^4+2*a*((phi1^2+theta1^2)*a^5+(2*rq*theta1^2-6*phi1^2-r1^2)*a^3+2*t1*phi1*(rq+4)*a^2+((-rq^2-8*rq)*phi1^2+rq^2*theta1^2-r1^2*rq-2*t1^2)*a+2*rq*t1*phi1*(rq+4))*r^3+4*(a*phi1-t1)*a*(a^4*phi1-(1/2)*a^3*t1+phi1*((7/2)*rq+1)*a^2-(1/2)*t1*(rq+2)*a+(5/2)*rq^2*phi1)*r^2-2*(a*phi1-t1)*(a^2+rq)*a*(phi1*(rq+2)*a^2-2*a*t1+rq^2*phi1)*r+a^2*(a*phi1-t1)^2*(a^2+rq)^2)*cos(theta)^2+2*a^2*sin(theta)*r^4*r1*theta1*(a^2+r^2-2*r+rq)*cos(theta)-((-phi1^2-theta1^2)*r^8+(4*phi1^2+4*theta1^2)*r^7-2*(phi1^2+theta1^2)*(a^2+rq+2)*r^6+((5*phi1^2+4*theta1^2)*a^2-2*a*t1*phi1+4*rq*phi1^2+4*rq*theta1^2-r1^2+t1^2)*r^5+((-phi1^2-theta1^2)*a^4+((-3*rq-4)*phi1^2-2*rq*theta1^2+r1^2)*a^2+2*t1*phi1*(rq+4)*a-rq^2*phi1^2-rq^2*theta1^2+(r1^2-t1^2)*rq-4*t1^2)*r^4+2*(a*phi1-t1)^2*(a^2+3*rq+2)*r^3-2*((rq+2)*a^2+rq^2+4*rq)*(a*phi1-t1)^2*r^2+(a*phi1-t1)^2*(a^2+rq)*(a^2+5*rq)*r-rq*(a*phi1-t1)^2*(a^2+rq)^2)*r)/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)^3), theta1, (sin(theta)*a^4*((phi1^2+theta1^2)*a^4+((2*phi1^2+theta1^2)*r^2+(-4*phi1^2-2*theta1^2)*r+2*rq*phi1^2+rq*theta1^2-r1^2)*a^2+phi1^2*(r^2-2*r+rq)^2)*cos(theta)^5-2*a^4*r*r1*theta1*(a^2+r^2-2*r+rq)*cos(theta)^4+2*sin(theta)*a^2*((phi1^2+theta1^2)*a^4+((2*phi1^2+theta1^2)*r^2+(-4*phi1^2-2*theta1^2)*r+2*rq*phi1^2+rq*theta1^2-r1^2)*a^2+phi1^2*(r^2-2*r+rq)^2)*r^2*cos(theta)^3-4*a^2*r^3*r1*theta1*(a^2+r^2-2*r+rq)*cos(theta)^2+2*sin(theta)*(phi1^2*(r-(1/2)*rq)*a^6-2*t1*(r-(1/2)*rq)*phi1*a^5+(((1/2)*theta1^2+(1/2)*phi1^2)*r^4+3*phi1^2*r^3-(3/2)*(rq+4/3)*phi1^2*r^2+(2*phi1^2*rq+t1^2)*r-(1/2)*rq*(phi1^2*rq+t1^2))*a^4-4*t1*(r-(1/2)*rq)*phi1*(r^2-r+(1/2)*rq)*a^3+((phi1^2+(1/2)*theta1^2)*r^6+(phi1^2-theta1^2)*r^5+((-(1/2)*rq-4)*phi1^2+(1/2)*rq*theta1^2-(1/2)*r1^2)*r^4+(4*phi1^2*rq+t1^2)*r^3+(-rq^2*phi1^2-(1/2)*t1^2*(rq+4))*r^2+2*r*rq*t1^2-(1/2)*rq^2*t1^2)*a^2-2*t1*(r-(1/2)*rq)*phi1*(r^2-2*r+rq)*r^2*a+(1/2)*phi1^2*r^6*(r^2-2*r+rq))*cos(theta)-2*r^5*r1*theta1*(a^2+r^2-2*r+rq))/((a^2+r^2-2*r+rq)*(r^2+a^2*cos(theta)^2)^3), phi1, (-2*a^4*phi1*theta1*(a^2+r^2-2*r+rq)*cos(theta)^5-2*a^4*sin(theta)*r1*phi1*(r-1)*cos(theta)^4-4*a^2*(a^2+r^2-2*r+rq)*phi1*(r^2-r+(1/2)*rq)*theta1*cos(theta)^3-2*sin(theta)*a^2*r1*(a^2*phi1-a*t1+2*phi1*(r-1/2)*r^2)*cos(theta)^2-4*(phi1*(r-(1/2)*rq)*a^2-t1*(r-(1/2)*rq)*a+(1/2)*phi1*r^4)*(a^2+r^2-2*r+rq)*theta1*cos(theta)+2*sin(theta)*r1*(phi1*(r-rq)*a^2-t1*(r-rq)*a-r^2*phi1*(r^2-2*r+rq))*r)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(a^2+r^2-2*r+rq))]';

endfunction

//For Kerr (Q=0<>a)
function Y=Kerr(V)//Q=0, J<>0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, (2*a^3*r1*phi1*(a-r)*(a+r)*cos(theta)^4+4*a^3*r*phi1*theta1*sin(theta)*(a^2+r^2-2*r)*cos(theta)^3-2*a*r1*(a^4*phi1+3*phi1*r^4-a^3*t1-a*r^2*t1)*cos(theta)^2-4*a^2*r*theta1*sin(theta)*(a^2+r^2-2*r)*(a*phi1-t1)*cos(theta)+2*r^2*r1*(a^3*phi1+3*a*phi1*r^2-a^2*t1-r^2*t1))/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2), r1, (-a^4*(phi1*r^2-2*phi1*r+a*(a*phi1-r1))*(phi1*r^2-2*phi1*r+a*(a*phi1+r1))*(r-1)*cos(theta)^6+2*a^6*r1*theta1*sin(theta)*(a^2+r^2-2*r)*cos(theta)^5+a^2*(-2*phi1^2*r^7+7*phi1^2*r^6+((-3*phi1^2+theta1^2)*a^2-4*phi1^2)*r^5+(-4*a^2*theta1^2+2*a*phi1*t1-4*phi1^2)*r^4+2*a*(a^3*theta1^2+(8*phi1^2+r1^2+2*theta1^2)*a-4*t1*phi1)*r^3+((-9*phi1^2-4*theta1^2)*a^4+4*a^3*t1*phi1+(-8*phi1^2-r1^2)*a^2+8*a*t1*phi1)*r^2+((phi1^2+theta1^2)*a^3+(8*phi1^2-r1^2)*a-8*t1*phi1)*a^3*r-2*phi1^2*a^6+2*t1*phi1*a^5)*cos(theta)^4+4*a^4*r^2*r1*theta1*sin(theta)*(a^2+r^2-2*r)*cos(theta)^3+(-phi1^2*r^9+4*phi1^2*r^8+(2*a^2*theta1^2-4*phi1^2)*r^7+((-2*phi1^2-8*theta1^2)*a^2-2*a*t1*phi1)*r^6+((3*phi1^2+4*theta1^2)*a^4+(r1^2+8*theta1^2)*a^2+8*a*t1*phi1)*r^5+((-3*phi1^2-8*theta1^2)*a^4-6*a^3*t1*phi1+(8*phi1^2+r1^2+t1^2)*a^2-8*a*t1*phi1)*r^4+2*((phi1^2+theta1^2)*a^4+(-6*phi1^2-r1^2)*a^2+8*a*t1*phi1-2*t1^2)*a^2*r^3+4*a^2*(a^3*phi1-(1/2)*t1*a^2+phi1*a-t1)*(a*phi1-t1)*r^2-4*a^4*(a*phi1-t1)^2*r+a^6*(a*phi1-t1)^2)*cos(theta)^2+2*a^2*r^4*r1*theta1*sin(theta)*(a^2+r^2-2*r)*cos(theta)-((-phi1^2-theta1^2)*r^7+(4*phi1^2+4*theta1^2)*r^6-2*(phi1^2+theta1^2)*(a^2+2)*r^5+((5*phi1^2+4*theta1^2)*a^2-2*a*t1*phi1-r1^2+t1^2)*r^4+((-phi1^2-theta1^2)*a^4+(-4*phi1^2+r1^2)*a^2+8*a*t1*phi1-4*t1^2)*r^3+2*(a^2+2)*(a*phi1-t1)^2*r^2-4*a^2*(a*phi1-t1)^2*r+a^4*(a*phi1-t1)^2)*r^2)/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^3), theta1, (a^4*(r^4*phi1^2-4*phi1^2*r^3+((2*phi1^2+theta1^2)*a^2+4*phi1^2)*r^2+(-4*phi1^2-2*theta1^2)*a^2*r+(phi1^2+theta1^2)*a^4-a^2*r1^2)*sin(theta)*cos(theta)^5-2*a^4*r*r1*theta1*(a^2+r^2-2*r)*cos(theta)^4+2*a^2*(r^4*phi1^2-4*phi1^2*r^3+((2*phi1^2+theta1^2)*a^2+4*phi1^2)*r^2+(-4*phi1^2-2*theta1^2)*a^2*r+(phi1^2+theta1^2)*a^4-a^2*r1^2)*sin(theta)*r^2*cos(theta)^3-4*a^2*r^3*r1*theta1*(a^2+r^2-2*r)*cos(theta)^2+2*((1/2)*phi1^2*r^7-phi1^2*r^6+a^2*(phi1^2+(1/2)*theta1^2)*r^5+((phi1^2-theta1^2)*a^2-2*a*t1*phi1)*r^4+(1/2)*((phi1^2+theta1^2)*a^3+(-8*phi1^2-r1^2)*a+8*t1*phi1)*a*r^3+(3*a^4*phi1^2-4*a^3*phi1*t1+a^2*t1^2)*r^2-2*a^2*(a*phi1-t1)^2*r+a^4*(a*phi1-t1)^2)*sin(theta)*r*cos(theta)-2*r^5*r1*theta1*(a^2+r^2-2*r))/((a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^3), phi1, (-2*a^4*phi1*theta1*(a^2+r^2-2*r)*cos(theta)^5-2*a^4*sin(theta)*r1*phi1*(r-1)*cos(theta)^4-4*a^2*r*phi1*theta1*(r-1)*(a^2+r^2-2*r)*cos(theta)^3-2*(2*phi1*r^3-phi1*r^2+a*(a*phi1-t1))*a^2*r1*sin(theta)*cos(theta)^2-4*theta1*(a^2+r^2-2*r)*((1/2)*phi1*r^3+a*(a*phi1-t1))*r*cos(theta)+2*(-phi1*r^3+2*phi1*r^2+a*(a*phi1-t1))*r1*sin(theta)*r^2)/(sin(theta)*(a^2+r^2-2*r)*(r^2+a^2*cos(theta)^2)^2)]';
endfunction

//For Reissner-Nordström (a=0<>Q)
function Y=ReissnerNordstrom(V)//Q<>0, J=0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, -2*(r-rq)*r1*t1/(r*(r^2-2*r+rq)), r1, (-r^4*phi1^2*(r^2-2*r+rq)^2*cos(theta)^2+(phi1^2+theta1^2)*r^8+(-4*phi1^2-4*theta1^2)*r^7+2*(phi1^2+theta1^2)*(rq+2)*r^6+((-4*phi1^2-4*theta1^2)*rq+r1^2-t1^2)*r^5+((phi1^2+theta1^2)*rq^2+(-r1^2+t1^2)*rq+4*t1^2)*r^4-6*(rq+2/3)*t1^2*r^3+2*rq*t1^2*(rq+4)*r^2-5*r*rq^2*t1^2+rq^3*t1^2)/(r^5*(r^2-2*r+rq)), theta1, (sin(theta)*cos(theta)*phi1^2*r-2*r1*theta1)/r, phi1, -2*r1*phi1/r-2*cos(theta)*theta1*phi1/sin(theta)]';
endfunction

//For Schwarzszchild (a=Q=0)
function Y=Schwarzschild(V)//Q=J=0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, -2*r1*t1/(r*(r-2)), r1, (-r^3*phi1^2*(r-2)^2*cos(theta)^2+(phi1^2+theta1^2)*r^5+(-4*phi1^2-4*theta1^2)*r^4+(4*phi1^2+4*theta1^2)*r^3+(r1^2-t1^2)*r^2+4*r*t1^2-4*t1^2)/(r^3*(r-2)), theta1, (sin(theta)*cos(theta)*phi1^2*r-2*r1*theta1)/r, phi1, -2*r1*phi1/r-2*cos(theta)*theta1*phi1/sin(theta)]';
endfunction

//For Minkowsky (a=Q=M=0)
function Y=Minkowsky(V)
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, 0, r1, r*(-cos(theta)^2*phi1^2+phi1^2+theta1^2), theta1, (sin(theta)*cos(theta)*phi1^2*r-2*theta1*r1)/r, phi1, -2*phi1*r1/r-2*cos(theta)*phi1*theta1/sin(theta)]';
endfunction








//The Cosmological (Lambda<>0) versions of the above functions (only those with a difference)

function g=cosmo_met_mat(V,rs,rq,a,Lambda)
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & a<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
    elseif (rq==0 & a<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
    elseif (rq<>0 & a==0) then
        chi=1;
        g=diag([(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)/(chi^2*r^2), -3*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), r^2, r^2*sin(theta)^2/chi^2]);
    else
        chi=1;
        g=diag([(1/3)*(Lambda*r^3-3*r+6)/(chi^2*r), -3*r/(Lambda*r^3-3*r+6), r^2, r^2*sin(theta)^2/chi^2]);
    end
endfunction


function [g,gi,drg,dthg]=cosmo_metric_matrix(V,rs,rq,a,Lambda)
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & a<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r+3*rq)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r+3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*sin(theta)^2)],4,4);
        drg=matrix([(1/3)*(2*Lambda*cos(theta)^4*a^4*r+4*Lambda*cos(theta)^2*a^2*r^3+2*Lambda*r^5+6*a^2*cos(theta)^2-6*r^2+6*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -(2/3)*sin(theta)^2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, (6*a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2+18*r*((1/3)*Lambda*r^4+a^2-r+rq))/(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)^2, 0, 0, 0, 0, 6*r/(Lambda*a^2*cos(theta)^2+3), 0, -(2/3)*sin(theta)^2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2/3)*sin(theta)^2*(a^4*(Lambda*a^2*r+3*r-3)*cos(theta)^4+((2*Lambda*r^3+3)*a^4+(6*r^3+3*r^2-3*r*rq)*a^2)*cos(theta)^2+r*((Lambda*r^4-3*r+3*rq)*a^2+3*r^4))/(chi^2*(r^2+a^2*cos(theta)^2)^2)],4,4);
        dthg=matrix([(2/3)*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -(2/3)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 6*a^2*cos(theta)*sin(theta)/(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq), 0, 0, 0, 0, 6*a^2*cos(theta)*sin(theta)*(Lambda*r^2-3)/(Lambda*a^2*cos(theta)^2+3)^2, 0, -(2/3)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^4*cos(theta)^4+2*r^2*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+(Lambda*r^4+6*r-3*rq)*a^4+r^2*(Lambda*r^4+3*r^2+12*r-6*rq)*a^2+3*r^6)*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^2)],4,4);
    elseif (rq==0 & a<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))/((r^2+a^2*cos(theta)^2)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*sin(theta)^2)],4,4);
        drg=matrix([(1/3)*(2*Lambda*cos(theta)^4*a^4*r+4*Lambda*cos(theta)^2*a^2*r^3+2*Lambda*r^5+6*a^2*cos(theta)^2-6*r^2)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -(2/3)*(Lambda*cos(theta)^4*a^4*r+a^2*(2*Lambda*r^3+3)*cos(theta)^2+Lambda*r^5-3*r^2)*a*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, (6*a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2+6*Lambda*r^5+18*a^2*r-18*r^2)/(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)^2, 0, 0, 0, 0, 6*r/(Lambda*a^2*cos(theta)^2+3), 0, -(2/3)*(Lambda*cos(theta)^4*a^4*r+a^2*(2*Lambda*r^3+3)*cos(theta)^2+Lambda*r^5-3*r^2)*a*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2/3)*(a^4*(Lambda*a^2*r+3*r-3)*cos(theta)^4+((2*Lambda*r^3+3)*a^4+(6*r^3+3*r^2)*a^2)*cos(theta)^2+(Lambda*r^5-3*r^2)*a^2+3*r^5)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)^2)],4,4);
        dthg=matrix([(2/3)*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, -(2/3)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 6*a^2*cos(theta)*sin(theta)/(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2), 0, 0, 0, 0, 6*a^2*cos(theta)*sin(theta)*(Lambda*r^2-3)/(Lambda*a^2*cos(theta)^2+3)^2, 0, -(2/3)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^2), 0, 0, (2/3)*cos(theta)*sin(theta)*(a^4*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^4+2*r^2*a^2*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^2+r*((Lambda*r^3+6)*a^4+r^2*(Lambda*r^3+3*r+12)*a^2+3*r^5))/(chi^2*(r^2+a^2*cos(theta)^2)^2)],4,4);
    elseif (rq<>0 & a==0) then
        chi=1;
        g=diag([(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)/(chi^2*r^2), -3*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), r^2, r^2*sin(theta)^2/chi^2]);
        gi=diag([3*chi^2*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), (1/3)*(-Lambda*r^4+3*r^2-6*r+3*rq)/r^2, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
        drg=diag([(1/3)*(2*Lambda*r^4-6*r+6*rq)/(chi^2*r^3), 6*r*(Lambda*r^4-3*r+3*rq)/(Lambda*r^4-3*r^2+6*r-3*rq)^2, 2*r, 2*r*sin(theta)^2/chi^2]);
        dthg=diag([0, 0, 0, 2*r^2*sin(theta)*cos(theta)/chi^2]);
    else
        chi=1;
        g=diag([(1/3)*(Lambda*r^3-3*r+6)/(chi^2*r), -3*r/(Lambda*r^3-3*r+6), r^2, r^2*sin(theta)^2/chi^2]);
        gi=diag([3*chi^2*r/(Lambda*r^3-3*r+6), (1/3)*(-Lambda*r^3+3*r-6)/r, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
        drg=diag([(1/3)*(2*Lambda*r^3-6)/(chi^2*r^2), 6*(Lambda*r^3-3)/(Lambda*r^3-3*r+6)^2, 2*r, 2*r*sin(theta)^2/chi^2]);
        dthg=diag([0, 0, 0, 2*r^2*sin(theta)*cos(theta)/chi^2]);
    end
endfunction

function gi=cosmo_inv_met_mat(V,rs,rq,a,Lambda)
    r=V(2); theta=V(3); phi=V(4);
    if (a<>0 & rq<>0) then
        chi=1+Lambda*a^2/3;
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r+3*rq)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r+3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*sin(theta)^2)],4,4);
    elseif (a<>0 & rq==0)
        chi=1+Lambda*a^2/3;
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))/((r^2+a^2*cos(theta)^2)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*sin(theta)^2)],4,4);
    elseif (a==0 & rq<>0)
        chi=1;
        gi=diag([3*chi^2*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), (1/3)*(-Lambda*r^4+3*r^2-6*r+3*rq)/r^2, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
    else
        chi=1;
        gi=diag([3*chi^2*r/(Lambda*r^3-3*r+6), (1/3)*(-Lambda*r^3+3*r-6)/r, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
    end
endfunction

function [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(V,rs,rq,a,Lambda)
    r=V(2); theta=V(3); phi=V(4);
    if (a<>0 & rq<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r+3*rq)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r+3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*sin(theta)^2)],4,4);
        drpg=matrix([-6*chi^2*(Lambda*a^4*(a^2+r^2)*((Lambda*a^2+3)*r^3-9*r^2+(Lambda*a^4+3*a^2+6*rq)*r+3*a^2)*cos(theta)^4+2*a^2*((Lambda^2*a^2+3*Lambda)*r^7-(15/2)*Lambda*r^6+2*Lambda*(Lambda*a^4+3*a^2+(9/4)*rq)*r^5+(-3*Lambda*a^2+9/2)*r^4+(Lambda^2*a^6+3*Lambda*a^4+3*Lambda*a^2*rq-18)*r^3+((9/2)*Lambda*a^4+9*a^2+18*rq)*r^2-(3/2)*rq*(Lambda*a^4+3*rq)*r+(9/2)*a^4)*cos(theta)^2+r*((Lambda^2*a^2+3*Lambda)*r^8+(2*Lambda^2*a^4+6*Lambda*a^2)*r^6+(12*Lambda*a^2-9)*r^5+(Lambda^2*a^6+3*Lambda*a^4-6*Lambda*a^2*rq+9*rq)*r^4+(12*Lambda*a^4-18*a^2)*r^3-6*a^2*(Lambda*a^2*rq-3*rq-6)*r^2+(-9*a^4-36*a^2*rq)*r+9*a^4*rq+9*a^2*rq^2))/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)^2), 0, 0, -6*chi^2*(Lambda*a^4*(Lambda*r^5+2*Lambda*a^2*r^3-3*r^2+(Lambda*a^4+3*rq)*r+3*a^2)*cos(theta)^4+2*(Lambda^2*r^7+2*Lambda^2*a^2*r^5+(3/2)*Lambda*r^4+Lambda^2*a^4*r^3+((9/2)*Lambda*a^2-9/2)*r^2-(3/2)*rq*(Lambda*a^2-3)*r+(9/2)*a^2)*a^2*cos(theta)^2+r*(r^8*Lambda^2+2*Lambda^2*a^2*r^6+12*Lambda*r^5+Lambda*(Lambda*a^4-6*rq)*r^4+(12*Lambda*a^2-27)*r^3+(-6*Lambda*a^2*rq+18*rq+36)*r^2+(-9*a^2-36*rq)*r+9*a^2*rq+9*rq^2))*a/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)^2), 0, (1/3)*(-2*a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2-6*r*((1/3)*Lambda*r^4+a^2-r+rq))/(r^2+a^2*cos(theta)^2)^2, 0, 0, 0, 0, -(2/3)*(Lambda*a^2*cos(theta)^2+3)*r/(r^2+a^2*cos(theta)^2)^2, 0, -6*chi^2*(Lambda*a^4*(Lambda*r^5+2*Lambda*a^2*r^3-3*r^2+(Lambda*a^4+3*rq)*r+3*a^2)*cos(theta)^4+2*(Lambda^2*r^7+2*Lambda^2*a^2*r^5+(3/2)*Lambda*r^4+Lambda^2*a^4*r^3+((9/2)*Lambda*a^2-9/2)*r^2-(3/2)*rq*(Lambda*a^2-3)*r+(9/2)*a^2)*a^2*cos(theta)^2+r*(r^8*Lambda^2+2*Lambda^2*a^2*r^6+12*Lambda*r^5+Lambda*(Lambda*a^4-6*rq)*r^4+(12*Lambda*a^2-27)*r^3+(-6*Lambda*a^2*rq+18*rq+36)*r^2+(-9*a^2-36*rq)*r+9*a^2*rq+9*rq^2))*a/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)^2), 0, 0, 6*(Lambda*a^6*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^6-a^4*(-3*Lambda^2*r^5-9*Lambda*r^2+(Lambda^2*a^4-3*Lambda*a^2+3*Lambda*rq+9)*r+3*Lambda*a^2-9)*cos(theta)^4-2*(((3/2)*a^2*Lambda^2-(9/2)*Lambda)*r^5+(Lambda^2*a^4-3*Lambda*a^2+9)*r^3+((9/2)*Lambda*a^2-27/2)*r^2-(3/2)*rq*(Lambda*a^2-3)*r+(9/2)*a^2)*a^2*cos(theta)^2-r*(r^8*Lambda^2+(2*Lambda^2*a^2-6*Lambda)*r^6+12*Lambda*r^5+(Lambda^2*a^4-3*Lambda*a^2-6*Lambda*rq+9)*r^4+(12*Lambda*a^2-36)*r^3+(-6*Lambda*a^2*rq+18*rq+36)*r^2+(-9*a^2-36*rq)*r+9*a^2*rq+9*rq^2))*chi^2/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)^2*sin(theta)^2)],4,4);
        dthpg=matrix([6*(Lambda*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^4*cos(theta)^4+2*Lambda*a^2*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*cos(theta)^2+Lambda^2*a^4*r^4+(Lambda*r^2+3)*(Lambda*r^4+6*r-3*rq)*a^2+3*r^2*(Lambda*r^4+6*r-3*rq))*chi^2*a^2*cos(theta)*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)^2*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 6*chi^2*(a^4*Lambda^2*(a^2+r^2)*cos(theta)^4+2*a^2*Lambda*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*cos(theta)^2+Lambda^2*a^2*r^4+Lambda^2*r^6+6*Lambda*r^3-3*Lambda*r^2*rq+18*r-9*rq)*a^3*cos(theta)*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, -(2/3)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*cos(theta)*a^2*sin(theta)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 0, 0, -(2/3)*sin(theta)*cos(theta)*a^2*(Lambda*r^2-3)/(r^2+a^2*cos(theta)^2)^2, 0, 6*chi^2*(a^4*Lambda^2*(a^2+r^2)*cos(theta)^4+2*a^2*Lambda*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*cos(theta)^2+Lambda^2*a^2*r^4+Lambda^2*r^6+6*Lambda*r^3-3*Lambda*r^2*rq+18*r-9*rq)*a^3*cos(theta)*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)), 0, 0, -6*chi^2*(Lambda^2*cos(theta)^8*a^8+(-2*Lambda^2*a^8+6*Lambda*a^6)*cos(theta)^6+a^4*(Lambda^2*a^4+(-3*Lambda^2*r^2-3*Lambda)*a^2-3*r^4*Lambda^2+9*Lambda*r^2-18*Lambda*r+9*rq*Lambda+9)*cos(theta)^4+2*(a^4*r^2*Lambda^2-3*Lambda*(r^2-2*r+rq)*a^2-(Lambda*r^2+3)*(Lambda*r^4-3*r^2+6*r-3*rq))*a^2*cos(theta)^2+Lambda^2*a^4*r^4+(Lambda^2*r^6-3*Lambda*r^4+6*Lambda*r^3-3*Lambda*r^2*rq+18*r-9*rq)*a^2-3*r^2*(Lambda*r^4-3*r^2+6*r-3*rq))*cos(theta)/((cos(theta)+1)*(r^2+a^2*cos(theta)^2)^2*(cos(theta)-1)*(Lambda*a^2*cos(theta)^2+3)^2*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)*sin(theta))],4,4);
    elseif (a<>0 & rq==0)
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        gi=matrix([3*chi^2*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))/((r^2+a^2*cos(theta)^2)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, (-Lambda*a^2*r^2-Lambda*r^4+3*a^2+3*r^2-6*r)/(3*a^2*cos(theta)^2+3*r^2), 0, 0, 0, 0, (Lambda*a^2*cos(theta)^2+3)/(3*a^2*cos(theta)^2+3*r^2), 0, 3*chi^2*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*a/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)), 0, 0, -3*chi^2*(a^4*cos(theta)^4*Lambda-Lambda*cos(theta)^2*a^4-Lambda*a^2*r^2-Lambda*r^4+3*a^2*cos(theta)^2+3*r^2-6*r)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*sin(theta)^2)],4,4);
        drpg=matrix([-6*chi^2*(Lambda*a^4*(a^2+r^2)*((Lambda*a^2+3)*r^3-9*r^2+(Lambda*a^4+3*a^2)*r+3*a^2)*cos(theta)^4+2*a^2*((Lambda^2*a^2+3*Lambda)*r^7-(15/2)*Lambda*r^6+(2*Lambda^2*a^4+6*Lambda*a^2)*r^5+(-3*Lambda*a^2+9/2)*r^4+(Lambda^2*a^6+3*Lambda*a^4-18)*r^3+((9/2)*Lambda*a^4+9*a^2)*r^2+(9/2)*a^4)*cos(theta)^2+((Lambda^2*a^2+3*Lambda)*r^7+(2*Lambda^2*a^4+6*Lambda*a^2)*r^5+(12*Lambda*a^2-9)*r^4+(Lambda^2*a^6+3*Lambda*a^4)*r^3+(12*Lambda*a^4-18*a^2)*r^2+36*a^2*r-9*a^4)*r^2)/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0, 0, -6*chi^2*a*(a^4*Lambda*(Lambda*a^4*r+2*Lambda*a^2*r^3+Lambda*r^5+3*a^2-3*r^2)*cos(theta)^4+(2*Lambda^2*a^2*r^7+4*Lambda^2*a^4*r^5+3*a^2*Lambda*r^4+2*Lambda^2*a^6*r^3+(9*Lambda*a^4-9*a^2)*r^2+9*a^4)*cos(theta)^2+r^2*(r^7*Lambda^2+2*Lambda^2*a^2*r^5+12*Lambda*r^4+a^4*r^3*Lambda^2+(12*Lambda*a^2-27)*r^2+36*r-9*a^2))/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0, (1/3)*(-2*a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2-2*Lambda*r^5-6*a^2*r+6*r^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 0, 0, -(2/3)*(Lambda*a^2*cos(theta)^2+3)*r/(r^2+a^2*cos(theta)^2)^2, 0, -6*chi^2*a*(a^4*Lambda*(Lambda*a^4*r+2*Lambda*a^2*r^3+Lambda*r^5+3*a^2-3*r^2)*cos(theta)^4+(2*Lambda^2*a^2*r^7+4*Lambda^2*a^4*r^5+3*a^2*Lambda*r^4+2*Lambda^2*a^6*r^3+(9*Lambda*a^4-9*a^2)*r^2+9*a^4)*cos(theta)^2+r^2*(r^7*Lambda^2+2*Lambda^2*a^2*r^5+12*Lambda*r^4+a^4*r^3*Lambda^2+(12*Lambda*a^2-27)*r^2+36*r-9*a^2))/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 6*chi^2*(Lambda*a^6*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^6-(-3*Lambda^2*r^5-9*Lambda*r^2+(Lambda^2*a^4-3*Lambda*a^2+9)*r+3*Lambda*a^2-9)*a^4*cos(theta)^4-2*(((3/2)*Lambda^2*a^2-(9/2)*Lambda)*r^5+(Lambda^2*a^4-3*Lambda*a^2+9)*r^3+((9/2)*Lambda*a^2-27/2)*r^2+(9/2)*a^2)*a^2*cos(theta)^2-r^2*(r^7*Lambda^2+(2*Lambda^2*a^2-6*Lambda)*r^5+12*Lambda*r^4+(Lambda^2*a^4-3*Lambda*a^2+9)*r^3+(12*Lambda*a^2-36)*r^2+36*r-9*a^2))/((r^2+a^2*cos(theta)^2)^2*sin(theta)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)^2*(Lambda*a^2*cos(theta)^2+3))],4,4);
        dthpg=matrix([6*chi^2*cos(theta)*(Lambda*a^4*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^4+2*r*Lambda*a^2*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3)*cos(theta)^2+r*(a^4*r^3*Lambda^2+(Lambda*r^2+3)*(Lambda*r^3+6)*a^2+3*Lambda*r^5+18*r^2))*a^2*sin(theta)/((r^2+a^2*cos(theta)^2)^2*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(Lambda*a^2*cos(theta)^2+3)^2), 0, 0, 6*chi^2*(a^4*Lambda^2*(a^2+r^2)*cos(theta)^4+2*a^2*r*Lambda*(Lambda*a^2*r+Lambda*r^3+6)*cos(theta)^2+Lambda^2*a^2*r^4+Lambda^2*r^6+6*Lambda*r^3+18*r)*cos(theta)*a^3*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)^2), 0, -(2/3)*sin(theta)*cos(theta)*a^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)/(r^2+a^2*cos(theta)^2)^2, 0, 0, 0, 0, -(2/3)*sin(theta)*cos(theta)*a^2*(Lambda*r^2-3)/(r^2+a^2*cos(theta)^2)^2, 0, 6*chi^2*(a^4*Lambda^2*(a^2+r^2)*cos(theta)^4+2*a^2*r*Lambda*(Lambda*a^2*r+Lambda*r^3+6)*cos(theta)^2+Lambda^2*a^2*r^4+Lambda^2*r^6+6*Lambda*r^3+18*r)*cos(theta)*a^3*sin(theta)/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*a^2*cos(theta)^2+3)^2), 0, 0, -6*chi^2*cos(theta)*(Lambda^2*cos(theta)^8*a^8+(-2*Lambda^2*a^8+6*Lambda*a^6)*cos(theta)^6+a^4*(Lambda^2*a^4+(-3*Lambda^2*r^2-3*Lambda)*a^2-3*r^4*Lambda^2+9*Lambda*r^2-18*r*Lambda+9)*cos(theta)^4+2*r*a^2*(a^4*r*Lambda^2-3*Lambda*(r-2)*a^2-Lambda^2*r^5-6*Lambda*r^2+9*r-18)*cos(theta)^2+r*(a^4*r^3*Lambda^2+(Lambda^2*r^5-3*Lambda*r^3+6*Lambda*r^2+18)*a^2-3*r^2*(Lambda*r^3-3*r+6)))/((r^2+a^2*cos(theta)^2)^2*(cos(theta)-1)*(cos(theta)+1)*sin(theta)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(Lambda*a^2*cos(theta)^2+3)^2)],4,4);
    elseif (a==0 & rq<>0)
        chi=1;
        g=diag([(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)/(chi^2*r^2), -3*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), r^2, r^2*sin(theta)^2/chi^2]);
        gi=diag([3*chi^2*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), (1/3)*(-Lambda*r^4+3*r^2-6*r+3*rq)/r^2, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
        drpg=diag([-6*chi^2*r*(Lambda*r^4-3*r+3*rq)/(Lambda*r^4-3*r^2+6*r-3*rq)^2, (1/3)*(-2*Lambda*r^4+6*r-6*rq)/r^3, -2/r^3, -2*chi^2/(r^3*sin(theta)^2)]);
        dthpg=diag([0, 0, 0, -2*chi^2*cos(theta)/(r^2*sin(theta)^3)]);
    else
        chi=1;
        g=diag([(1/3)*(Lambda*r^3-3*r+6)/(chi^2*r), -3*r/(Lambda*r^3-3*r+6), r^2, r^2*sin(theta)^2/chi^2]);
        gi=diag([3*chi^2*r/(Lambda*r^3-3*r+6), (1/3)*(-Lambda*r^3+3*r-6)/r, 1/r^2, chi^2/(r^2*sin(theta)^2)]);
        drpg=diag([-6*chi^2*(Lambda*r^3-3)/(Lambda*r^3-3*r+6)^2, (1/3)*(-2*Lambda*r^3+6)/r^2, -2/r^3, -2*chi^2/(r^3*sin(theta)^2)]);
        dthpg=diag([0, 0, 0, -2*chi^2*cos(theta)/(r^2*sin(theta)^3)]);
    end
endfunction

function [g,Gam]=cosmo_metric_with_christoffel(V)
    G=1; c=1;
    r=V(2); theta=V(3); phi=V(4);
    if (rq<>0 & k<>0 & rs<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        Gam=list(matrix([0, 2*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)*(a^2+r^2)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)), -2*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)^2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6*a*(a^2*(a^2-r^2+r*rq)*cos(theta)^2-a^2*r^2+a^2*r*rq-3*r^4+2*r^3*rq)*sin(theta)^2/((Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*(r^2+a^2*cos(theta)^2)^2), (12*r-6*rq)*a^3*sin(theta)^3*cos(theta)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)^2), 0],4,4)', matrix([(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, -(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*sin(theta)^2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, (-a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2-3*r*((1/3)*Lambda*r^4+a^2-r+rq))/((Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*(r^2+a^2*cos(theta)^2)), -2*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2), 0, 0, 0, r*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3)), 0, -(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*sin(theta)^2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, (1/9)*((-3+(Lambda*a^2+3)*r)*a^4*cos(theta)^4+((2*Lambda*a^4+6*a^2)*r^3+3*a^2*r^2-3*a^2*r*rq+3*a^4)*cos(theta)^2+r*((Lambda*a^2+3)*r^4-3*a^2*r+3*a^2*rq))*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)^3)],4,4)', matrix([-(1/9)*(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, (1/9)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)*(Lambda*a^2*cos(theta)^2+3)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, -(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)/(((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)*(r^2+a^2*cos(theta)^2)), 0, 0, 0, 2*r/(r^2+a^2*cos(theta)^2), sin(theta)*cos(theta)*a^2*(Lambda*r^2-3)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)), 0, (1/9)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)*(Lambda*a^2*cos(theta)^2+3)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, -(1/9)*(Lambda*a^2*cos(theta)^2+3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^4*cos(theta)^4+2*r^2*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^2*cos(theta)^2+(Lambda*r^4+6*r-3*rq)*a^4+r^2*(Lambda*r^4+3*r^2+12*r-6*rq)*a^2+3*r^6)*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^3)],4,4)', matrix([0, 2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2+3*r*rq)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)), -2*cos(theta)*a*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2*a^4*(Lambda*r^3-3*r+3)*cos(theta)^4-6*(-(2/3)*Lambda*r^5+2*r^3+a^2-r^2)*a^2*cos(theta)^2+2*Lambda*r^7-6*r^5+6*a^2*r^2-6*a^2*r*rq+12*r^4-6*r^3*rq)/((Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*(r^2+a^2*cos(theta)^2)^2), 2*cos(theta)*(Lambda*cos(theta)^6*a^6+a^4*(2*Lambda*r^2+3)*cos(theta)^4+a^2*(Lambda*r^4+6*r^2-6*r+3*rq)*cos(theta)^2+(6*r-3*rq)*a^2+3*r^4)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0],4,4)');
    elseif (rq==0 & k<>0 & rs<>0) then
        chi=1+Lambda*a^2/3;
        g=matrix([(1/3)*(-a^4*cos(theta)^4*Lambda+Lambda*cos(theta)^2*a^4+Lambda*a^2*r^2+Lambda*r^4-3*a^2*cos(theta)^2-3*r^2+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, (-3*a^2*cos(theta)^2-3*r^2)/(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r), 0, 0, 0, 0, (3*a^2*cos(theta)^2+3*r^2)/(Lambda*a^2*cos(theta)^2+3), 0, -(1/3)*a*(a^2*Lambda*(a^2+r^2)*cos(theta)^2+Lambda*a^2*r^2+Lambda*r^4+6*r)*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2)), 0, 0, (1/3)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*a^2*cos(theta)^2+r*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3))*sin(theta)^2/(chi^2*(r^2+a^2*cos(theta)^2))],4,4);
        Gam=list(matrix([0, 2*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2)*(a^2+r^2)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)), -2*sin(theta)*cos(theta)*a^2*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)^2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6*sin(theta)^2*a*((a^4-a^2*r^2)*cos(theta)^2-a^2*r^2-3*r^4)/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)), 12*sin(theta)^3*cos(theta)*a^3*r/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)^2), 0],4,4)', matrix([(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, -(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*a*sin(theta)^2*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, (-a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2-Lambda*r^5-3*a^2*r+3*r^2)/((r^2+a^2*cos(theta)^2)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)), -2*a^2*cos(theta)*sin(theta)/(r^2+a^2*cos(theta)^2), 0, 0, 0, r*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)), 0, -(1/9)*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)*a*sin(theta)^2*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, (1/9)*((-3+(Lambda*a^2+3)*r)*a^4*cos(theta)^4+((2*Lambda*a^4+6*a^2)*r^3+3*a^2*r^2+3*a^4)*cos(theta)^2+(Lambda*a^2+3)*r^5-3*a^2*r^2)*sin(theta)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)/(chi^2*(r^2+a^2*cos(theta)^2)^3)],4,4)', matrix([-(1/9)*(Lambda*a^2*cos(theta)^2+3)*sin(theta)*cos(theta)*a^2*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, (1/9)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)*(Lambda*a^2*cos(theta)^2+3)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, -(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)/((r^2+a^2*cos(theta)^2)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)), 0, 0, 0, 2*r/(r^2+a^2*cos(theta)^2), sin(theta)*cos(theta)*a^2*(Lambda*r^2-3)/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2)), 0, (1/9)*(a^2+r^2)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)*(Lambda*a^2*cos(theta)^2+3)*a*cos(theta)*sin(theta)/(chi^2*(r^2+a^2*cos(theta)^2)^3), 0, 0, -(1/9)*cos(theta)*(a^4*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^4+2*r^2*a^2*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^2+r*((Lambda*r^3+6)*a^4+r^2*(Lambda*r^3+3*r+12)*a^2+3*r^5))*sin(theta)*(Lambda*a^2*cos(theta)^2+3)/(chi^2*(r^2+a^2*cos(theta)^2)^3)],4,4)', matrix([0, 2*a*(Lambda*cos(theta)^4*a^4*r+2*Lambda*cos(theta)^2*a^2*r^3+Lambda*r^5+3*a^2*cos(theta)^2-3*r^2)/((r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r)), -2*cos(theta)*a*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (2*a^4*(Lambda*r^3-3*r+3)*cos(theta)^4-6*(-(2/3)*Lambda*r^5+2*r^3+a^2-r^2)*a^2*cos(theta)^2+2*Lambda*r^7-6*r^5+6*a^2*r^2+12*r^4)/((r^2+a^2*cos(theta)^2)^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)), 2*cos(theta)*(Lambda*cos(theta)^6*a^6+a^4*(2*Lambda*r^2+3)*cos(theta)^4+a^2*r*(Lambda*r^3+6*r-6)*cos(theta)^2+3*r^4+6*a^2*r)/(sin(theta)*(r^2+a^2*cos(theta)^2)^2*(Lambda*a^2*cos(theta)^2+3)), 0],4,4)');
    elseif (rq<>0 & k==0 & rs<>0) then
        chi=1;
        g=diag([(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)/(chi^2*r^2), -3*r^2/(Lambda*r^4-3*r^2+6*r-3*rq), r^2, r^2*sin(theta)^2/chi^2]);
        Gam=list(matrix([0, (2*Lambda*r^4-6*r+6*rq)/(r*(Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],4,4)', matrix([(1/9)*(Lambda*r^4-3*r^2+6*r-3*rq)*(Lambda*r^4-3*r+3*rq)/(r^5*chi^2), 0, 0, 0, 0, (-Lambda*r^4+3*r-3*rq)/(r*(Lambda*r^4-3*r^2+6*r-3*rq)), 0, 0, 0, 0, (1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)/r, 0, 0, 0, 0, (1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)*sin(theta)^2/(r*chi^2)],4,4)', matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 2/r, 0, 0, 0, 0, 0, -sin(theta)*cos(theta)/chi^2],4,4)', matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2/r, 2*cos(theta)/sin(theta), 0],4,4)');
    elseif (rq==0 & k==0 & rs<>0) then
        chi=1;
        g=diag([(1/3)*(Lambda*r^3-3*r+6)/(chi^2*r), -3*r/(Lambda*r^3-3*r+6), r^2, r^2*sin(theta)^2/chi^2]);
        Gam=list(matrix([0, (2*Lambda*r^3-6)/((Lambda*r^3-3*r+6)*r), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],4,4)', matrix([(1/9)*(Lambda*r^3-3*r+6)*(Lambda*r^3-3)/(r^3*chi^2), 0, 0, 0, 0, (-Lambda*r^3+3)/((Lambda*r^3-3*r+6)*r), 0, 0, 0, 0, (1/3)*Lambda*r^3-r+2, 0, 0, 0, 0, (1/3)*(Lambda*r^3-3*r+6)*sin(theta)^2/chi^2],4,4)', matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 2/r, 0, 0, 0, 0, 0, -sin(theta)*cos(theta)/chi^2],4,4)', matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2/r, 2*cos(theta)/sin(theta), 0],4,4)');
    elseif rs==0 then
        g=diag([-1,1,1,1]); Gam=list(zeros(4,4),zeros(4,4),zeros(4,4),zeros(4,4));
    end
endfunction











//Hamilton equations

function X=cosmo_init_conds_hamiltonian(V)
    g=cosmo_met_mat(V,rs,rq,a,Lambda);
    X=[V(1:4)',(g*V(5:8))']';
endfunction

function Y=cosmo_Hamilton_equations(V)
    [g,gi,drpg,dthpg]=cosmo_inverse_metric_matrix(V(1:4),rs,rq,a,Lambda);
    Y=[gi(1,:)*V(5:8),gi(2,:)*V(5:8),gi(3,:)*V(5:8),gi(4,:)*V(5:8),0,-sum(V(5:8).*(drpg*V(5:8)))/2,-sum(V(5:8).*(dthpg*V(5:8)))/2,0]';
endfunction










//Carter-Newman (inspired by PYYY)

function Y=cosmo_Carter_Newman(V)
    ep=1e-10; chi=1+Lambda*a^2/3;
    r=V(1); th=V(2); ph=V(3); pr=V(4); pth=V(5);
    Dr=(1-Lambda*r^2/3)*(r^2+a^2)-2*r+rq; Dt=1+Lambda*a^2*cos(th)^2/3; S=r^2+a^2*cos(th)^2;
    Drp=-2/3*Lambda*a^2*r-4/3*Lambda*r^3+2*r-2; Dtp=-2/3*Lambda*a^2*cos(th)*sin(th);
    if abs(Dr)<ep then
        Dr=sign(Dr)*ep;
    end
    if abs(S)<ep then
        S=sign(S)*ep;
    end
    if abs(sin(th))<ep then
        th=sign(th)*asin(ep);
    end
    Pofr=chi*(E*(r^2+a^2)-a*Lz)+e*r*sqrt(rq); Wofth=chi*(a*E*sin(th)-Lz/sin(th));
    prp=(((2*chi*E*r+e*sqrt(rq))*Pofr-Drp*(k-mu*r^2)/2)/Dr+mu*r-Drp*pr^2)/S;
    pthp=((Dtp*(k+mu*a^2*cos(th)^2)/2-chi^2*cos(th)*sin(th)*(a^2*E^2-Lz^2/sin(th)^4))/Dt-mu*a^2*cos(th)*sin(th)-Dtp*pth^2)/S;
    rp=Dr*pr/S;
    thp=Dt*pth/S;
    php=chi/S*(a*Pofr/Dr-Wofth/(Dt*sin(th)));
    Y=[rp,thp,php,prp,pthp]';
endfunction














//Euler-Lagrange equation (computed with simple Christoffel symbols in Maple)


function Y=cosmo_KerrNewman(V)//Q<>0, J<>0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    chi=1+Lambda*a^2/3;
    Y=[t1, (-2*sin(theta)*t1*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)*a^8*Lambda*theta1*cos(theta)^7-2*r1*a^6*Lambda*((a^3*r*t1*Lambda+3*phi1*a^2+r^3*t1*Lambda*a-3*phi1*r*(r-rq))*a*sin(theta)^2-r*t1*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq))*cos(theta)^6-6*sin(theta)*t1*a^6*Lambda*theta1*((Lambda*r^4-3*r^2+2*r-rq)*a^2+r^2*(Lambda*r^4-3*r^2+8*r-4*rq))*cos(theta)^5-6*r1*a^4*((t1*Lambda*(Lambda*r^3+1)*a^3+phi1*(Lambda*r*rq+3)*a^2+r*t1*Lambda*(Lambda*r^4+3*r-rq)*a-4*r*(Lambda*r^3-(3/4)*Lambda*r^2*rq+(3/4)*r-(3/4)*rq)*phi1)*a*sin(theta)^2-t1*((Lambda^2*r^3+Lambda)*a^4+(Lambda^2*r^5+3*Lambda*r^3+3*Lambda*r^2-Lambda*r*rq+3)*a^2+3*Lambda*r^5-4*Lambda*r^4+2*r^3*rq*Lambda+3*r^2-6*r+3*rq))*cos(theta)^4-12*sin(theta)*a^4*theta1*((Lambda*t1*a^3+(-Lambda*phi1*r^2+3*phi1)*a^2+a*r^2*t1*Lambda-phi1*(Lambda*r^4-3*r^2+6*r-3*rq))*a*(r-(1/2)*rq)*sin(theta)^2-(Lambda*(r-(1/2)*rq)*a^4-(1/2)*(Lambda*r^2-3)*(Lambda*r^4+2*r-rq)*a^2-(1/2)*r^8*Lambda^2+(3/2)*Lambda*r^6-5*Lambda*r^5+(5/2)*r^4*rq*Lambda+3*r^3+(-(3/2)*rq-12)*r^2+12*r*rq-3*rq^2)*t1)*cos(theta)^3-6*r1*a^2*((r*t1*Lambda*(Lambda*r^4+rq)*a^3-phi1*r*(Lambda*r^3-Lambda*r^2*rq-3*rq)*a^2+t1*(Lambda^2*r^7+4*Lambda*r^4-Lambda*r^3*rq+6*r-3*rq)*a-3*r^3*phi1*(Lambda*r^3-(2/3)*Lambda*r^2*rq+4*r-3*rq))*a*sin(theta)^2-t1*(r*Lambda*(Lambda*r^4+rq)*a^4+(Lambda^2*r^7+3*Lambda*r^5+4*Lambda*r^4-r^3*rq*Lambda+(3*rq+6)*r-3*rq)*a^2+3*Lambda*r^7-2*Lambda*r^6+Lambda*r^5*rq+(3*rq+6)*r^3-9*r^2*rq+3*r*rq^2))*cos(theta)^2-12*sin(theta)*((a^3*r^2*t1*Lambda-r^2*phi1*(Lambda*r^2-3)*a^2+t1*(Lambda*r^4+6*r-3*rq)*a-r^2*phi1*(Lambda*r^4-3*r^2+6*r-3*rq))*a*(r-(1/2)*rq)*sin(theta)^2-t1*(r^2*Lambda*(r-(1/2)*rq)*a^4-(1/6)*(Lambda*r^4+6*r-3*rq)*(Lambda*r^4-3*r^2-6*r+3*rq)*a^2-(1/6)*r^2*(Lambda*r^4+6*r-3*rq)*(Lambda*r^4-3*r^2+6*r-3*rq)))*a^2*theta1*cos(theta)-2*r1*r*(a*(r^2*t1*Lambda*(Lambda*r^4-3*r+3*rq)*a^3-9*r^2*phi1*(r-rq)*a^2+t1*(Lambda*r^4-3*r+3*rq)*(Lambda*r^4+6*r-3*rq)*a-27*r^4*(r-(2/3)*rq)*phi1)*sin(theta)^2-t1*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)*(Lambda*r^4-3*r+3*rq)))/((r^2+a^2*cos(theta)^2)*(a^6*Lambda*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*cos(theta)^6+(a^2*Lambda^2*(a^2+r^2)^2*sin(theta)^2-a^6*Lambda^2+(Lambda^2*r^4+3*Lambda*r^2+12*Lambda*r-6*Lambda*rq+9)*a^2+3*Lambda*r^4+9*r^2-18*r+9*rq)*a^4*cos(theta)^4+2*(Lambda*a^2*(a^2+r^2)*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2-a^6*r^2*Lambda^2-(3/2)*Lambda*(Lambda*r^4+4*r-2*rq)*a^4+(-(1/2)*Lambda^2*r^6-(3/2)*Lambda*r^4+9*r^2)*a^2-(3/2)*Lambda*r^6+3*Lambda*r^5+(-(3/2)*rq*Lambda+9)*r^4-18*r^3+(9*rq+18)*r^2-18*r*rq+(9/2)*rq^2)*a^2*cos(theta)^2+a^2*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)^2*sin(theta)^2-(Lambda*a^2*r^2+Lambda*r^4-3*r^2+6*r-3*rq)*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4))), r1, -(1/9)*(Lambda*cos(theta)^4*a^4*r+(2*Lambda*a^2*r^3+3*a^2)*cos(theta)^2+r*(Lambda*r^4-3*r+3*rq))*t1^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)/((r^2+a^2*cos(theta)^2)^3*chi^2)+(2/9)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*sin(theta)^2*(Lambda*cos(theta)^4*a^4*r+(2*Lambda*a^2*r^3+3*a^2)*cos(theta)^2+r*(Lambda*r^4-3*r+3*rq))*a*phi1*t1/((r^2+a^2*cos(theta)^2)^3*chi^2)+r1^2*(a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2+3*r*((1/3)*Lambda*r^4+a^2-r+rq))/((Lambda*a^2*r^2+Lambda*r^4-3*a^2-3*r^2+6*r-3*rq)*(r^2+a^2*cos(theta)^2))+2*a^2*cos(theta)*sin(theta)*theta1*r1/(r^2+a^2*cos(theta)^2)-theta1^2*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2-3*rq)*r/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3))-(1/9)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq)*(a^4*(Lambda*a^2*r+3*r-3)*cos(theta)^4+2*a^2*((Lambda*r^3+3/2)*a^2+3*r*(r^2+(1/2)*r-(1/2)*rq))*cos(theta)^2+r*((Lambda*r^4-3*r+3*rq)*a^2+3*r^4))*phi1^2*sin(theta)^2/((r^2+a^2*cos(theta)^2)^3*chi^2), theta1, (1/9)*(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r-3*rq)*t1^2/((r^2+a^2*cos(theta)^2)^3*chi^2)-(2/9)*(Lambda*a^2*cos(theta)^2+3)*(a^4*Lambda*(a^2+r^2)*cos(theta)^4+2*a^2*(Lambda*a^2*r^2+Lambda*r^4+3*r-(3/2)*rq)*cos(theta)^2+6*(r-(1/2)*rq)*a^2*sin(theta)^2+r^2*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq))*sin(theta)*a*phi1*cos(theta)*t1/((r^2+a^2*cos(theta)^2)^3*chi^2)+(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)*r1^2/((r^2+a^2*cos(theta)^2)*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r-3*rq))-2*r*r1*theta1/(r^2+a^2*cos(theta)^2)-a^2*cos(theta)*sin(theta)*(Lambda*r^2-3)*theta1^2/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2))+(1/9)*cos(theta)*phi1^2*(Lambda*a^2*cos(theta)^2+3)*sin(theta)*((Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*a^4*cos(theta)^4+2*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+3*r-(3/2)*rq)*a^2+3*r^2*(r^2-r+(1/2)*rq))*a^2*cos(theta)^2+6*(a^2+r^2)*(r-(1/2)*rq)*a^2*sin(theta)^2+r^2*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4))/((r^2+a^2*cos(theta)^2)^3*chi^2), phi1, (-2*a^8*Lambda*theta1*(a^4*phi1*Lambda-Lambda*t1*a^3+phi1*(Lambda*r^2+3)*a^2-a*r^2*t1*Lambda+3*phi1*(r^2-2*r+rq))*cos(theta)^9-2*r1*sin(theta)*(a^2*r*phi1*Lambda-r*t1*Lambda*a+3*phi1*(r-1))*a^8*Lambda*cos(theta)^8-2*a^6*theta1*(a^2*Lambda^2*(a^2+r^2)*(a^2*phi1+phi1*r^2-a*t1)*sin(theta)^2-a^6*phi1*Lambda^2+a^5*t1*Lambda^2+a^4*r^2*phi1*Lambda^2-t1*Lambda*(Lambda*r^2+3)*a^3+2*(r^4*Lambda^2+3*r^2*Lambda+6*r*Lambda-3*rq*Lambda+9/2)*phi1*a^2-2*t1*(Lambda*r^4+(3/2)*r^2+3*r-(3/2)*rq)*Lambda*a+6*phi1*(Lambda*r^4-Lambda*r^3+((1/2)*rq*Lambda+3/2)*r^2-3*r+(3/2)*rq))*cos(theta)^7-2*r1*sin(theta)*(a^2*r*phi1*Lambda^2*(a^2+r^2)*sin(theta)^2-a^4*r*phi1*Lambda^2+2*phi1*Lambda*(Lambda*r^3+3)*a^2-3*t1*Lambda*(Lambda*r^3+r+1)*a+6*(-3/2+Lambda*r^3+(1/2)*r^2*Lambda+(-(1/2)*rq*Lambda+3/2)*r)*phi1)*a^6*cos(theta)^6-6*(a^2*Lambda*(a^4*r^2*phi1*Lambda-a^3*r^2*t1*Lambda+2*(Lambda*r^4+3*r-(3/2)*rq)*phi1*a^2-t1*(Lambda*r^4+4*r-2*rq)*a+r^2*phi1*(Lambda*r^4+6*r-3*rq))*sin(theta)^2-a^6*r^2*phi1*Lambda^2+r^2*t1*Lambda^2*a^5-phi1*Lambda*(Lambda*r^4+4*r-2*rq)*a^4+t1*Lambda*(Lambda*r^4-3*r^2+4*r-2*rq)*a^3+4*r^2*(r*Lambda-(1/2)*rq*Lambda+9/4)*phi1*a^2-3*t1*(Lambda*r^4+2*r-rq)*a+2*phi1*(Lambda*r^5+(-(1/2)*rq*Lambda+9/2)*r^4-9*r^3+((9/2)*rq+6)*r^2-6*r*rq+(3/2)*rq^2))*a^4*theta1*cos(theta)^5-6*r1*sin(theta)*(((Lambda*r^3+1)*a^2+r*(Lambda*r^4+3*r-rq))*a^2*Lambda*phi1*sin(theta)^2-phi1*Lambda*(Lambda*r^3+1)*a^4-3*phi1*(Lambda*r^2-Lambda*r*rq-1)*a^2-t1*(Lambda^2*r^5+3*Lambda*r^3+Lambda*r*rq+3)*a+phi1*(Lambda*r^4+9*r^3-6*r^2+6*r-3*rq))*a^4*cos(theta)^4-12*a^2*theta1*(a^4*(a^2+r^2)*Lambda*(r-(1/2)*rq)*phi1*sin(theta)^4+(1/2)*(phi1*Lambda*(Lambda*r^4-2*r+rq)*a^4-a^3*r^4*t1*Lambda^2+2*(Lambda^2*r^6+3*Lambda*r^3-(3/2)*Lambda*r^2*rq+3*r-(3/2)*rq)*phi1*a^2-t1*(Lambda^2*r^6+6*Lambda*r^3-3*Lambda*r^2*rq+6*r-3*rq)*a+(r^8*Lambda^2+8*Lambda*r^5-4*r^4*rq*Lambda+6*r^3+(-3*rq+12)*r^2-12*r*rq+3*rq^2)*phi1)*a^2*sin(theta)^2-(1/2)*a^6*r^4*phi1*Lambda^2+(1/2)*a^5*r^4*t1*Lambda^2-(5/6)*r^2*Lambda*phi1*(Lambda*r^4+(24/5)*r-(12/5)*rq)*a^4+(5/6)*r^2*t1*Lambda*(Lambda*r^4-(9/5)*r^2+(24/5)*r-(12/5)*rq)*a^3-(1/3)*(r^8*Lambda^2+3*Lambda*r^6+6*Lambda*r^5+(-3*rq*Lambda-27/2)*r^4+18*r^2-18*r*rq+(9/2)*rq^2)*phi1*a^2+(1/3)*t1*(r^8*Lambda^2-(9/2)*Lambda*r^6+9*Lambda*r^5-(9/2)*r^4*rq*Lambda-18*r^3+(9*rq+18)*r^2-18*r*rq+(9/2)*rq^2)*a-r^2*phi1*(Lambda*r^6-Lambda*r^5+((1/2)*rq*Lambda-9/2)*r^4+9*r^3+(-(9/2)*rq-6)*r^2+6*r*rq-(3/2)*rq^2))*cos(theta)^3-6*r1*sin(theta)*(a^2*(r*Lambda*(Lambda*r^4+rq)*a^2+Lambda^2*r^7+4*Lambda*r^4-r^3*rq*Lambda+6*r-3*rq)*phi1*sin(theta)^2-r*phi1*Lambda*(Lambda*r^4+rq)*a^4-(2/3)*(Lambda^2*r^7+9*Lambda*r^4-(9/2)*r^3*rq*Lambda+(9-(9/2)*rq)*r-(9/2)*rq)*phi1*a^2-(1/3)*r*t1*(Lambda^2*r^6+9*Lambda*r^4-3*Lambda*r^3+3*Lambda*r^2*rq+9*rq)*a-2*r*phi1*(Lambda*r^6+(1/2)*Lambda*r^5+(-(1/2)*rq*Lambda-9/2)*r^4+(9/2)*r^3+(-(3/2)*rq+3)*r^2-(9/2)*r*rq+(3/2)*rq^2))*a^2*cos(theta)^2-12*theta1*(a^4*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*(r-(1/2)*rq)*phi1*sin(theta)^4+(1/6)*(r^2*phi1*Lambda*(Lambda*r^4-6*r+3*rq)*a^4-a^3*r^6*t1*Lambda^2+2*(r^8*Lambda^2+9*r^3+(-(9/2)*rq-18)*r^2+18*r*rq-(9/2)*rq^2)*phi1*a^2-r^2*t1*(Lambda^2*r^6+6*Lambda*r^3-3*Lambda*r^2*rq+18*r-9*rq)*a+r^4*phi1*(Lambda^2*r^6+6*Lambda*r^3-3*Lambda*r^2*rq+18*r-9*rq))*a^2*sin(theta)^2-(1/6)*r^2*(a^4*r^2*phi1*Lambda-a^3*r^2*t1*Lambda+phi1*(Lambda*r^4+3*r^2+6*r-3*rq)*a^2-t1*(Lambda*r^4+6*r-3*rq)*a+3*r^4*phi1)*(Lambda*a^2*r^2+Lambda*r^4-3*r^2+6*r-3*rq))*cos(theta)-2*r1*sin(theta)*r*(a^2*phi1*(Lambda*r^4-3*r+3*rq)*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2-r^2*phi1*Lambda*(Lambda*r^4-3*r+3*rq)*a^4-(r^8*Lambda^2+3*Lambda*r^5+9*r^3+(-9*rq-18)*r^2+27*r*rq-9*rq^2)*phi1*a^2-3*r^2*t1*(Lambda*r^4-3*r+3*rq)*a-3*r^4*phi1*(Lambda*r^4-3*r^2+6*r-3*rq)))/((r^2+a^2*cos(theta)^2)*sin(theta)*(a^6*Lambda*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r+3*rq)*cos(theta)^6+(a^2*Lambda^2*(a^2+r^2)^2*sin(theta)^2-a^6*Lambda^2+(Lambda^2*r^4+3*Lambda*r^2+12*Lambda*r-6*Lambda*rq+9)*a^2+3*Lambda*r^4+9*r^2-18*r+9*rq)*a^4*cos(theta)^4+2*(Lambda*a^2*(a^2+r^2)*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)*sin(theta)^2-a^6*r^2*Lambda^2-(3/2)*Lambda*(Lambda*r^4+4*r-2*rq)*a^4+(-(1/2)*Lambda^2*r^6-(3/2)*Lambda*r^4+9*r^2)*a^2-(3/2)*Lambda*r^6+3*Lambda*r^5+(-(3/2)*rq*Lambda+9)*r^4-18*r^3+(9*rq+18)*r^2-18*r*rq+(9/2)*rq^2)*a^2*cos(theta)^2+a^2*(Lambda*a^2*r^2+Lambda*r^4+6*r-3*rq)^2*sin(theta)^2-(Lambda*a^2*r^2+Lambda*r^4-3*r^2+6*r-3*rq)*(a^4*Lambda*r^2+(Lambda*r^4+3*r^2+6*r-3*rq)*a^2+3*r^4)))]';
endfunction

function Y=cosmo_Kerr(V)//Q=0, J<>0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    chi=1+Lambda*a^2/3;
    Y=[t1, (-2*sin(theta)*t1*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*a^8*theta1*Lambda*cos(theta)^7-2*r1*a^6*(a*(Lambda*a^3*r*t1+Lambda*a*r^3*t1+3*a^2*phi1-3*phi1*r^2)*sin(theta)^2-r*t1*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r))*Lambda*cos(theta)^6-6*r*sin(theta)*((Lambda*r^3-3*r+2)*a^2+r^2*(Lambda*r^3-3*r+8))*t1*a^6*theta1*Lambda*cos(theta)^5-6*r1*(a*(t1*Lambda*(Lambda*r^3+1)*a^3+3*phi1*a^2+Lambda*r^2*t1*(Lambda*r^3+3)*a-4*r^4*phi1*Lambda-3*r^2*phi1)*sin(theta)^2-t1*((Lambda^2*r^3+Lambda)*a^4+(Lambda^2*r^5+3*Lambda*r^3+3*Lambda*r^2+3)*a^2+3*Lambda*r^5-4*Lambda*r^4+3*r^2-6*r))*a^4*cos(theta)^4-12*(a*(Lambda*t1*a^3+(-Lambda*phi1*r^2+3*phi1)*a^2+a*r^2*t1*Lambda-r*phi1*(Lambda*r^3-3*r+6))*sin(theta)^2-t1*(Lambda*a^4-(1/2)*(Lambda*r^3+2)*(Lambda*r^2-3)*a^2-(1/2)*r*(Lambda^2*r^6-3*Lambda*r^4+10*Lambda*r^3-6*r+24)))*r*sin(theta)*a^4*theta1*cos(theta)^3-6*r1*r*((a^3*r^4*t1*Lambda^2-r^3*phi1*Lambda*a^2+t1*(Lambda^2*r^6+4*Lambda*r^3+6)*a-3*r^3*phi1*(Lambda*r^2+4))*a*sin(theta)^2-t1*(Lambda^2*a^4*r^4+(Lambda^2*r^6+3*Lambda*r^4+4*Lambda*r^3+6)*a^2+3*Lambda*r^6-2*Lambda*r^5+6*r^2))*a^2*cos(theta)^2-12*r^2*sin(theta)*a^2*theta1*((a^3*r*t1*Lambda-r*phi1*(Lambda*r^2-3)*a^2+t1*(Lambda*r^3+6)*a-r^2*phi1*(Lambda*r^3-3*r+6))*a*sin(theta)^2-t1*(a^4*Lambda*r+(-(1/6)*Lambda^2*r^6+(1/2)*Lambda*r^4+3*r+6)*a^2-(1/6)*r^2*(Lambda*r^3+6)*(Lambda*r^3-3*r+6)))*cos(theta)-2*r1*r^3*(a*(r*t1*Lambda*(Lambda*r^3-3)*a^3-9*a^2*r*phi1+t1*(Lambda*r^3+6)*(Lambda*r^3-3)*a-27*r^3*phi1)*sin(theta)^2-t1*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3)*(Lambda*r^3-3)))/((a^6*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*Lambda*cos(theta)^6+a^4*(a^2*Lambda^2*(a^2+r^2)^2*sin(theta)^2-a^6*Lambda^2+(Lambda^2*r^4+3*Lambda*r^2+12*Lambda*r+9)*a^2+3*Lambda*r^4+9*r^2-18*r)*cos(theta)^4+2*(Lambda*a^2*(a^2+r^2)*(Lambda*a^2*r+Lambda*r^3+6)*sin(theta)^2-a^6*r*Lambda^2+(-(3/2)*r^3*Lambda^2-6*Lambda)*a^4+(-(1/2)*r^5*Lambda^2-(3/2)*Lambda*r^3+9*r)*a^2-(3/2)*Lambda*r^5+3*Lambda*r^4+9*r^3-18*r^2+18*r)*r*a^2*cos(theta)^2+r^2*(a^2*(Lambda*a^2*r+Lambda*r^3+6)^2*sin(theta)^2-(Lambda*a^2*r+Lambda*r^3-3*r+6)*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3)))*(r^2+a^2*cos(theta)^2)), r1, -(1/9)*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(Lambda*cos(theta)^4*a^4*r+(2*Lambda*a^2*r^3+3*a^2)*cos(theta)^2+Lambda*r^5-3*r^2)*t1^2/((r^2+a^2*cos(theta)^2)^3*chi^2)+(2/9)*sin(theta)^2*a*t1*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*phi1*(Lambda*cos(theta)^4*a^4*r+(2*Lambda*a^2*r^3+3*a^2)*cos(theta)^2+Lambda*r^5-3*r^2)/((r^2+a^2*cos(theta)^2)^3*chi^2)+(a^2*(3+2*Lambda*r^3+(Lambda*a^2-3)*r)*cos(theta)^2+Lambda*r^5+3*a^2*r-3*r^2)*r1^2/((Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*(r^2+a^2*cos(theta)^2))+2*a^2*cos(theta)*sin(theta)*theta1*r1/(r^2+a^2*cos(theta)^2)-r*(Lambda*r^4+(Lambda*a^2-3)*r^2+6*r-3*a^2)*theta1^2/((r^2+a^2*cos(theta)^2)*(Lambda*a^2*cos(theta)^2+3))-(1/9)*sin(theta)^2*phi1^2*(a^4*(Lambda*a^2*r+3*r-3)*cos(theta)^4+((2*Lambda*r^3+3)*a^4+(6*r^3+3*r^2)*a^2)*cos(theta)^2+r^2*((Lambda*r^3-3)*a^2+3*r^3))*((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)/((r^2+a^2*cos(theta)^2)^3*chi^2), theta1, (1/9)*(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)*(a^4*cos(theta)^4*Lambda+2*Lambda*cos(theta)^2*a^2*r^2+Lambda*r^4+6*r)*t1^2/((r^2+a^2*cos(theta)^2)^3*chi^2)-(2/9)*sin(theta)*a*t1*cos(theta)*(a^4*Lambda*(a^2+r^2)*cos(theta)^4+2*a^2*r*(Lambda*a^2*r+Lambda*r^3+3)*cos(theta)^2+a^2*Lambda*r^4+Lambda*r^6+6*sin(theta)^2*a^2*r+6*r^3)*(Lambda*a^2*cos(theta)^2+3)*phi1/((r^2+a^2*cos(theta)^2)^3*chi^2)+(Lambda*a^2*cos(theta)^2+3)*a^2*cos(theta)*sin(theta)*r1^2/(((Lambda*r^2-3)*a^2+Lambda*r^4-3*r^2+6*r)*(r^2+a^2*cos(theta)^2))-2*r*r1*theta1/(r^2+a^2*cos(theta)^2)-a^2*cos(theta)*sin(theta)*(Lambda*r^2-3)*theta1^2/((Lambda*a^2*cos(theta)^2+3)*(r^2+a^2*cos(theta)^2))+(1/9)*sin(theta)*(Lambda*a^2*cos(theta)^2+3)*phi1^2*cos(theta)*(a^4*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*cos(theta)^4+2*a^2*r*(a^4*Lambda*r+(Lambda*r^3+3*r+3)*a^2+3*r^3-3*r^2)*cos(theta)^2+((6*a^4+6*a^2*r^2)*sin(theta)^2+(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3)*r^2)*r)/((r^2+a^2*cos(theta)^2)^3*chi^2), phi1, (-2*a^8*theta1*(a^4*phi1*Lambda-Lambda*t1*a^3+phi1*(Lambda*r^2+3)*a^2-a*r^2*t1*Lambda+3*r*phi1*(r-2))*Lambda*cos(theta)^9-2*r1*sin(theta)*(a^2*r*phi1*Lambda-r*t1*Lambda*a+3*phi1*(r-1))*a^8*Lambda*cos(theta)^8-2*a^6*theta1*(a^2*Lambda^2*(a^2+r^2)*(a^2*phi1+phi1*r^2-a*t1)*sin(theta)^2-a^6*phi1*Lambda^2+a^5*t1*Lambda^2+a^4*r^2*phi1*Lambda^2-t1*Lambda*(Lambda*r^2+3)*a^3+2*phi1*(r^4*Lambda^2+3*r^2*Lambda+6*r*Lambda+9/2)*a^2-2*r*t1*(Lambda*r^3+(3/2)*r+3)*Lambda*a+6*r*phi1*(Lambda*r^3-r^2*Lambda+(3/2)*r-3))*cos(theta)^7-2*r1*sin(theta)*a^6*(a^2*r*phi1*Lambda^2*(a^2+r^2)*sin(theta)^2-a^4*r*phi1*Lambda^2+2*phi1*Lambda*(Lambda*r^3+3)*a^2-3*t1*Lambda*(Lambda*r^3+r+1)*a+6*(Lambda*r^3+(1/2)*r^2*Lambda+(3/2)*r-3/2)*phi1)*cos(theta)^6-6*r*(a^2*(a^4*r*phi1*Lambda-a^3*r*t1*Lambda+(2*Lambda*phi1*r^3+6*phi1)*a^2+(-Lambda*r^3*t1-4*t1)*a+r^2*phi1*(Lambda*r^3+6))*Lambda*sin(theta)^2-a^6*r*phi1*Lambda^2+r*t1*Lambda^2*a^5-phi1*Lambda*(Lambda*r^3+4)*a^4+t1*Lambda*(Lambda*r^3-3*r+4)*a^3+(4*Lambda*phi1*r^2+9*phi1*r)*a^2+(-3*Lambda*r^3*t1-6*t1)*a+2*r*phi1*(Lambda*r^3+(9/2)*r^2-9*r+6))*a^4*theta1*cos(theta)^5-6*r1*sin(theta)*a^4*(((Lambda*r^3+1)*a^2+Lambda*r^5+3*r^2)*phi1*a^2*Lambda*sin(theta)^2-phi1*Lambda*(Lambda*r^3+1)*a^4+(-3*Lambda*phi1*r^2+3*phi1)*a^2-t1*(Lambda^2*r^5+3*Lambda*r^3+3)*a+r*phi1*(Lambda*r^3+9*r^2-6*r+6))*cos(theta)^4-12*r*a^2*(a^4*phi1*Lambda*(a^2+r^2)*sin(theta)^4+(1/2)*a^2*(phi1*Lambda*(Lambda*r^3-2)*a^4-a^3*r^3*t1*Lambda^2+2*phi1*(Lambda^2*r^5+3*Lambda*r^2+3)*a^2-t1*(Lambda^2*r^5+6*Lambda*r^2+6)*a+r*phi1*(Lambda^2*r^6+8*Lambda*r^3+6*r+12))*sin(theta)^2-(1/2)*r*(a^6*r^2*phi1*Lambda^2-r^2*t1*Lambda^2*a^5+((5/3)*r^4*phi1*Lambda^2+8*r*phi1*Lambda)*a^4-(5/3)*r*t1*(Lambda*r^3-(9/5)*r+24/5)*Lambda*a^3+(2/3)*(Lambda^2*r^6+3*Lambda*r^4+6*Lambda*r^3-(27/2)*r^2+18)*phi1*a^2-(2/3)*t1*(Lambda^2*r^6-(9/2)*Lambda*r^4+9*Lambda*r^3-18*r+18)*a+2*r^2*phi1*(Lambda*r^4-Lambda*r^3-(9/2)*r^2+9*r-6)))*theta1*cos(theta)^3-6*r1*r*sin(theta)*a^2*(a^2*phi1*(Lambda^2*a^2*r^4+Lambda^2*r^6+4*Lambda*r^3+6)*sin(theta)^2-a^4*r^4*phi1*Lambda^2-(2/3)*phi1*(Lambda^2*r^6+9*Lambda*r^3+9)*a^2-(1/3)*r^3*t1*Lambda*(Lambda*r^3+9*r-3)*a-2*r^2*phi1*(Lambda*r^4+(1/2)*Lambda*r^3-(9/2)*r^2+(9/2)*r+3))*cos(theta)^2-12*r^2*(a^4*phi1*(Lambda*a^2*r+Lambda*r^3+6)*sin(theta)^4+(1/6)*a^2*(r*phi1*Lambda*(Lambda*r^3-6)*a^4-a^3*r^4*t1*Lambda^2+2*phi1*(Lambda^2*r^6+9*r-18)*a^2-t1*r*(Lambda^2*r^5+6*Lambda*r^2+18)*a+r^3*phi1*(Lambda^2*r^5+6*Lambda*r^2+18))*sin(theta)^2-(1/6)*(Lambda*a^2*r+Lambda*r^3-3*r+6)*r^2*(a^4*r*phi1*Lambda-a^3*r*t1*Lambda+phi1*(Lambda*r^3+3*r+6)*a^2+(-Lambda*r^3*t1-6*t1)*a+3*r^3*phi1))*theta1*cos(theta)-2*r1*r^3*sin(theta)*(a^2*phi1*(Lambda*r^3-3)*(Lambda*a^2*r+Lambda*r^3+6)*sin(theta)^2-r*phi1*Lambda*(Lambda*r^3-3)*a^4-phi1*(Lambda^2*r^6+3*Lambda*r^3+9*r-18)*a^2-3*t1*r*(Lambda*r^3-3)*a-3*r^3*phi1*(Lambda*r^3-3*r+6)))/(sin(theta)*(a^6*(Lambda*a^4+(Lambda*r^2+3)*a^2+3*r^2-6*r)*Lambda*cos(theta)^6+a^4*(a^2*Lambda^2*(a^2+r^2)^2*sin(theta)^2-a^6*Lambda^2+(Lambda^2*r^4+3*Lambda*r^2+12*Lambda*r+9)*a^2+3*Lambda*r^4+9*r^2-18*r)*cos(theta)^4+2*(Lambda*a^2*(a^2+r^2)*(Lambda*a^2*r+Lambda*r^3+6)*sin(theta)^2-a^6*r*Lambda^2+(-(3/2)*r^3*Lambda^2-6*Lambda)*a^4+(-(1/2)*r^5*Lambda^2-(3/2)*Lambda*r^3+9*r)*a^2-(3/2)*Lambda*r^5+3*Lambda*r^4+9*r^3-18*r^2+18*r)*r*a^2*cos(theta)^2+r^2*(a^2*(Lambda*a^2*r+Lambda*r^3+6)^2*sin(theta)^2-(Lambda*a^2*r+Lambda*r^3-3*r+6)*(a^4*Lambda*r+(Lambda*r^3+3*r+6)*a^2+3*r^3)))*(r^2+a^2*cos(theta)^2))]';
endfunction

function Y=cosmo_ReissnerNordstrom(V)//Q<>0, J=0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    chi=1;
    Y=[t1, -2*(Lambda*r^4-3*r+3*rq)*r1*t1/(r*(Lambda*r^4-3*r^2+6*r-3*rq)), r1, (1/9)*(-Lambda*r^4+3*r^2-6*r+3*rq)*(Lambda*r^4-3*r+3*rq)*t1^2/(r^5*chi^2)-9*((1/3)*Lambda*r^4+rq-r)*(-(1/3)*Lambda*r^4+r^2+rq-2*r)*r1^2/(r*(-Lambda*r^4+3*r^2-6*r+3*rq)^2)-(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)*theta1^2/r-(1/3)*(Lambda*r^4-3*r^2+6*r-3*rq)*sin(theta)^2*phi1^2/(r*chi^2), theta1, -2*r1*theta1/r+sin(theta)*cos(theta)*phi1^2/chi^2, phi1, -2*r1*phi1/r-2*cos(theta)*theta1*phi1/sin(theta)]';
endfunction

function Y=cosmo_Schwarzschild(V)//Q=J=0
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    chi=1;
    Y=[t1, -2*(Lambda*r^3-3)*r1*t1/(r*(Lambda*r^3-3*r+6)), r1, (1/9)*(-3*phi1^2*r^3*(Lambda*r^3-3*r+6)^2*sin(theta)^2+(-3*Lambda^2*chi^2*theta1^2-Lambda^3*t1^2)*r^9+(18*Lambda*chi^2*theta1^2+6*Lambda^2*t1^2)*r^7+(-36*Lambda*chi^2*theta1^2-9*Lambda^2*t1^2)*r^6+((9*chi^2*r1^2-9*t1^2)*Lambda-27*chi^2*theta1^2)*r^5+(108*chi^2*theta1^2+18*Lambda*t1^2)*r^4-108*theta1^2*chi^2*r^3+(-27*chi^2*r1^2+27*t1^2)*r^2-108*r*t1^2+108*t1^2)/(r^3*(Lambda*r^3-3*r+6)*chi^2), theta1, -2*r1*theta1/r+sin(theta)*cos(theta)*phi1^2/chi^2, phi1, -2*r1*phi1/r-2*cos(theta)*theta1*phi1/sin(theta)]';
endfunction

function Y=deSitter(V)
    t=V(1); t1=V(2); r=V(3); r1=V(4); theta=V(5); theta1=V(6); phi=V(7); phi1=V(8);
    Y=[t1, 0, r1, r*(-cos(theta)^2*phi1^2+phi1^2+theta1^2), theta1, (sin(theta)*cos(theta)*phi1^2*r-2*theta1*r1)/r, phi1, -2*phi1*r1/r-2*cos(theta)*phi1*theta1/sin(theta)]';
endfunction













//RGB values for the temperature of a black-body radiation (see http://www.vendian.org/mncharity/dir3/blackbody/ )
 blackbody=[ 1000  1.0000 0.0337 0.0000;
  1000  1.0000 0.0401 0.0000;
  1100  1.0000 0.0592 0.0000;
  1100  1.0000 0.0631 0.0000;
  1200  1.0000 0.0846 0.0000;
  1200  1.0000 0.0860 0.0000;
  1300  1.0000 0.1096 0.0000;
  1300  1.0000 0.1085 0.0000;
  1400  1.0000 0.1341 0.0000;
  1400  1.0000 0.1303 0.0000;
  1500  1.0000 0.1578 0.0000;
  1500  1.0000 0.1515 0.0000;
  1600  1.0000 0.1806 0.0000;
  1600  1.0000 0.1718 0.0000;
  1700  1.0000 0.2025 0.0000;
  1700  1.0000 0.1912 0.0000;
  1800  1.0000 0.2235 0.0000;
  1800  1.0000 0.2097 0.0000;
  1900  1.0000 0.2434 0.0000;
  1900  1.0000 0.2272 0.0000;
  2000  1.0000 0.2647 0.0033;
  2000  1.0000 0.2484 0.0061;
  2100  1.0000 0.2889 0.0120;
  2100  1.0000 0.2709 0.0153;
  2200  1.0000 0.3126 0.0219;
  2200  1.0000 0.2930 0.0257;
  2300  1.0000 0.3360 0.0331;
  2300  1.0000 0.3149 0.0373;
  2400  1.0000 0.3589 0.0454;
  2400  1.0000 0.3364 0.0501;
  2500  1.0000 0.3814 0.0588;
  2500  1.0000 0.3577 0.0640;
  2600  1.0000 0.4034 0.0734;
  2600  1.0000 0.3786 0.0790;
  2700  1.0000 0.4250 0.0889;
  2700  1.0000 0.3992 0.0950;
  2800  1.0000 0.4461 0.1054;
  2800  1.0000 0.4195 0.1119;
  2900  1.0000 0.4668 0.1229;
  2900  1.0000 0.4394 0.1297;
  3000  1.0000 0.4870 0.1411;
  3000  1.0000 0.4589 0.1483;
  3100  1.0000 0.5067 0.1602;
  3100  1.0000 0.4781 0.1677;
  3200  1.0000 0.5259 0.1800;
  3200  1.0000 0.4970 0.1879;
  3300  1.0000 0.5447 0.2005;
  3300  1.0000 0.5155 0.2087;
  3400  1.0000 0.5630 0.2216;
  3400  1.0000 0.5336 0.2301;
  3500  1.0000 0.5809 0.2433;
  3500  1.0000 0.5515 0.2520;
  3600  1.0000 0.5983 0.2655;
  3600  1.0000 0.5689 0.2745;
  3700  1.0000 0.6153 0.2881;
  3700  1.0000 0.5860 0.2974;
  3800  1.0000 0.6318 0.3112;
  3800  1.0000 0.6028 0.3207;
  3900  1.0000 0.6480 0.3346;
  3900  1.0000 0.6193 0.3444;
  4000  1.0000 0.6636 0.3583;
  4000  1.0000 0.6354 0.3684;
  4100  1.0000 0.6789 0.3823;
  4100  1.0000 0.6511 0.3927;
  4200  1.0000 0.6938 0.4066;
  4200  1.0000 0.6666 0.4172;
  4300  1.0000 0.7083 0.4310;
  4300  1.0000 0.6817 0.4419;
  4400  1.0000 0.7223 0.4556;
  4400  1.0000 0.6966 0.4668;
  4500  1.0000 0.7360 0.4803;
  4500  1.0000 0.7111 0.4919;
  4600  1.0000 0.7494 0.5051;
  4600  1.0000 0.7253 0.5170;
  4700  1.0000 0.7623 0.5299;
  4700  1.0000 0.7392 0.5422;
  4800  1.0000 0.7750 0.5548;
  4800  1.0000 0.7528 0.5675;
  4900  1.0000 0.7872 0.5797;
  4900  1.0000 0.7661 0.5928;
  5000  1.0000 0.7992 0.6045;
  5000  1.0000 0.7792 0.6180;
  5100  1.0000 0.8108 0.6293;
  5100  1.0000 0.7919 0.6433;
  5200  1.0000 0.8221 0.6541;
  5200  1.0000 0.8044 0.6685;
  5300  1.0000 0.8330 0.6787;
  5300  1.0000 0.8167 0.6937;
  5400  1.0000 0.8437 0.7032;
  5400  1.0000 0.8286 0.7187;
  5500  1.0000 0.8541 0.7277;
  5500  1.0000 0.8403 0.7437;
  5600  1.0000 0.8642 0.7519;
  5600  1.0000 0.8518 0.7686;
  5700  1.0000 0.8740 0.7760;
  5700  1.0000 0.8630 0.7933;
  5800  1.0000 0.8836 0.8000;
  5800  1.0000 0.8740 0.8179;
  5900  1.0000 0.8929 0.8238;
  5900  1.0000 0.8847 0.8424;
  6000  1.0000 0.9019 0.8473;
  6000  1.0000 0.8952 0.8666;
  6100  1.0000 0.9107 0.8707;
  6100  1.0000 0.9055 0.8907;
  6200  1.0000 0.9193 0.8939;
  6200  1.0000 0.9156 0.9147;
  6300  1.0000 0.9276 0.9168;
  6300  1.0000 0.9254 0.9384;
  6400  1.0000 0.9357 0.9396;
  6400  1.0000 0.9351 0.9619;
  6500  1.0000 0.9436 0.9621;
  6500  1.0000 0.9445 0.9853;
  6600  1.0000 0.9513 0.9844;
  6600  0.9917 0.9458 1.0000;
  6700  0.9937 0.9526 1.0000;
  6700  0.9696 0.9336 1.0000;
  6800  0.9726 0.9395 1.0000;
  6800  0.9488 0.9219 1.0000;
  6900  0.9526 0.9270 1.0000;
  6900  0.9290 0.9107 1.0000;
  7000  0.9337 0.9150 1.0000;
  7000  0.9102 0.9000 1.0000;
  7100  0.9157 0.9035 1.0000;
  7100  0.8923 0.8897 1.0000;
  7200  0.8986 0.8925 1.0000;
  7200  0.8753 0.8799 1.0000;
  7300  0.8823 0.8819 1.0000;
  7300  0.8591 0.8704 1.0000;
  7400  0.8668 0.8718 1.0000;
  7400  0.8437 0.8614 1.0000;
  7500  0.8520 0.8621 1.0000;
  7500  0.8289 0.8527 1.0000;
  7600  0.8379 0.8527 1.0000;
  7600  0.8149 0.8443 1.0000;
  7700  0.8244 0.8437 1.0000;
  7700  0.8014 0.8363 1.0000;
  7800  0.8115 0.8351 1.0000;
  7800  0.7885 0.8285 1.0000;
  7900  0.7992 0.8268 1.0000;
  7900  0.7762 0.8211 1.0000;
  8000  0.7874 0.8187 1.0000;
  8000  0.7644 0.8139 1.0000;
  8100  0.7761 0.8110 1.0000;
  8100  0.7531 0.8069 1.0000;
  8200  0.7652 0.8035 1.0000;
  8200  0.7423 0.8002 1.0000;
  8300  0.7548 0.7963 1.0000;
  8300  0.7319 0.7938 1.0000;
  8400  0.7449 0.7894 1.0000;
  8400  0.7219 0.7875 1.0000;
  8500  0.7353 0.7827 1.0000;
  8500  0.7123 0.7815 1.0000;
  8600  0.7260 0.7762 1.0000;
  8600  0.7030 0.7757 1.0000;
  8700  0.7172 0.7699 1.0000;
  8700  0.6941 0.7700 1.0000;
  8800  0.7086 0.7638 1.0000;
  8800  0.6856 0.7645 1.0000;
  8900  0.7004 0.7579 1.0000;
  8900  0.6773 0.7593 1.0000;
  9000  0.6925 0.7522 1.0000;
  9000  0.6693 0.7541 1.0000;
  9100  0.6848 0.7467 1.0000;
  9100  0.6617 0.7492 1.0000;
  9200  0.6774 0.7414 1.0000;
  9200  0.6543 0.7444 1.0000;
  9300  0.6703 0.7362 1.0000;
  9300  0.6471 0.7397 1.0000;
  9400  0.6635 0.7311 1.0000;
  9400  0.6402 0.7352 1.0000;
  9500  0.6568 0.7263 1.0000;
  9500  0.6335 0.7308 1.0000;
  9600  0.6504 0.7215 1.0000;
  9600  0.6271 0.7265 1.0000;
  9700  0.6442 0.7169 1.0000;
  9700  0.6208 0.7224 1.0000;
  9800  0.6382 0.7124 1.0000;
  9800  0.6148 0.7183 1.0000;
  9900  0.6324 0.7081 1.0000;
  9900  0.6089 0.7144 1.0000;
 10000  0.6268 0.7039 1.0000;
 10000  0.6033 0.7106 1.0000;
 10100  0.6213 0.6998 1.0000;
 10100  0.5978 0.7069 1.0000;
 10200  0.6161 0.6958 1.0000;
 10200  0.5925 0.7033 1.0000;
 10300  0.6109 0.6919 1.0000;
 10300  0.5873 0.6998 1.0000;
 10400  0.6060 0.6881 1.0000;
 10400  0.5823 0.6964 1.0000;
 10500  0.6012 0.6844 1.0000;
 10500  0.5774 0.6930 1.0000;
 10600  0.5965 0.6808 1.0000;
 10600  0.5727 0.6898 1.0000;
 10700  0.5919 0.6773 1.0000;
 10700  0.5681 0.6866 1.0000;
 10800  0.5875 0.6739 1.0000;
 10800  0.5637 0.6836 1.0000;
 10900  0.5833 0.6706 1.0000;
 10900  0.5593 0.6806 1.0000;
 11000  0.5791 0.6674 1.0000;
 11000  0.5551 0.6776 1.0000;
 11100  0.5750 0.6642 1.0000;
 11100  0.5510 0.6748 1.0000;
 11200  0.5711 0.6611 1.0000;
 11200  0.5470 0.6720 1.0000;
 11300  0.5673 0.6581 1.0000;
 11300  0.5432 0.6693 1.0000;
 11400  0.5636 0.6552 1.0000;
 11400  0.5394 0.6666 1.0000;
 11500  0.5599 0.6523 1.0000;
 11500  0.5357 0.6640 1.0000;
 11600  0.5564 0.6495 1.0000;
 11600  0.5322 0.6615 1.0000;
 11700  0.5530 0.6468 1.0000;
 11700  0.5287 0.6590 1.0000;
 11800  0.5496 0.6441 1.0000;
 11800  0.5253 0.6566 1.0000;
 11900  0.5463 0.6415 1.0000;
 11900  0.5220 0.6542 1.0000;
 12000  0.5431 0.6389 1.0000;
 12000  0.5187 0.6519 1.0000;
 12100  0.5400 0.6364 1.0000;
 12100  0.5156 0.6497 1.0000;
 12200  0.5370 0.6340 1.0000;
 12200  0.5125 0.6474 1.0000;
 12300  0.5340 0.6316 1.0000;
 12300  0.5095 0.6453 1.0000;
 12400  0.5312 0.6293 1.0000;
 12400  0.5066 0.6432 1.0000;
 12500  0.5283 0.6270 1.0000;
 12500  0.5037 0.6411 1.0000;
 12600  0.5256 0.6247 1.0000;
 12600  0.5009 0.6391 1.0000;
 12700  0.5229 0.6225 1.0000;
 12700  0.4982 0.6371 1.0000;
 12800  0.5203 0.6204 1.0000;
 12800  0.4955 0.6351 1.0000;
 12900  0.5177 0.6183 1.0000;
 12900  0.4929 0.6332 1.0000;
 13000  0.5152 0.6162 1.0000;
 13000  0.4904 0.6314 1.0000;
 13100  0.5128 0.6142 1.0000;
 13100  0.4879 0.6295 1.0000;
 13200  0.5104 0.6122 1.0000;
 13200  0.4854 0.6277 1.0000;
 13300  0.5080 0.6103 1.0000;
 13300  0.4831 0.6260 1.0000;
 13400  0.5057 0.6084 1.0000;
 13400  0.4807 0.6243 1.0000;
 13500  0.5035 0.6065 1.0000;
 13500  0.4785 0.6226 1.0000;
 13600  0.5013 0.6047 1.0000;
 13600  0.4762 0.6209 1.0000;
 13700  0.4991 0.6029 1.0000;
 13700  0.4740 0.6193 1.0000;
 13800  0.4970 0.6012 1.0000;
 13800  0.4719 0.6177 1.0000;
 13900  0.4950 0.5994 1.0000;
 13900  0.4698 0.6161 1.0000;
 14000  0.4930 0.5978 1.0000;
 14000  0.4677 0.6146 1.0000;
 14100  0.4910 0.5961 1.0000;
 14100  0.4657 0.6131 1.0000;
 14200  0.4891 0.5945 1.0000;
 14200  0.4638 0.6116 1.0000;
 14300  0.4872 0.5929 1.0000;
 14300  0.4618 0.6102 1.0000;
 14400  0.4853 0.5913 1.0000;
 14400  0.4599 0.6087 1.0000;
 14500  0.4835 0.5898 1.0000;
 14500  0.4581 0.6073 1.0000;
 14600  0.4817 0.5882 1.0000;
 14600  0.4563 0.6060 1.0000;
 14700  0.4799 0.5868 1.0000;
 14700  0.4545 0.6046 1.0000;
 14800  0.4782 0.5853 1.0000;
 14800  0.4527 0.6033 1.0000;
 14900  0.4765 0.5839 1.0000;
 14900  0.4510 0.6020 1.0000;
 15000  0.4749 0.5824 1.0000;
 15000  0.4493 0.6007 1.0000;
 15100  0.4733 0.5811 1.0000;
 15100  0.4477 0.5994 1.0000;
 15200  0.4717 0.5797 1.0000;
 15200  0.4460 0.5982 1.0000;
 15300  0.4701 0.5784 1.0000;
 15300  0.4445 0.5970 1.0000;
 15400  0.4686 0.5770 1.0000;
 15400  0.4429 0.5958 1.0000;
 15500  0.4671 0.5757 1.0000;
 15500  0.4413 0.5946 1.0000;
 15600  0.4656 0.5745 1.0000;
 15600  0.4398 0.5935 1.0000;
 15700  0.4641 0.5732 1.0000;
 15700  0.4384 0.5923 1.0000;
 15800  0.4627 0.5720 1.0000;
 15800  0.4369 0.5912 1.0000;
 15900  0.4613 0.5708 1.0000;
 15900  0.4355 0.5901 1.0000;
 16000  0.4599 0.5696 1.0000;
 16000  0.4341 0.5890 1.0000;
 16100  0.4586 0.5684 1.0000;
 16100  0.4327 0.5879 1.0000;
 16200  0.4572 0.5673 1.0000;
 16200  0.4313 0.5869 1.0000;
 16300  0.4559 0.5661 1.0000;
 16300  0.4300 0.5859 1.0000;
 16400  0.4546 0.5650 1.0000;
 16400  0.4287 0.5848 1.0000;
 16500  0.4534 0.5639 1.0000;
 16500  0.4274 0.5838 1.0000;
 16600  0.4521 0.5628 1.0000;
 16600  0.4261 0.5829 1.0000;
 16700  0.4509 0.5617 1.0000;
 16700  0.4249 0.5819 1.0000;
 16800  0.4497 0.5607 1.0000;
 16800  0.4236 0.5809 1.0000;
 16900  0.4485 0.5597 1.0000;
 16900  0.4224 0.5800 1.0000;
 17000  0.4474 0.5586 1.0000;
 17000  0.4212 0.5791 1.0000;
 17100  0.4462 0.5576 1.0000;
 17100  0.4201 0.5781 1.0000;
 17200  0.4451 0.5566 1.0000;
 17200  0.4189 0.5772 1.0000;
 17300  0.4440 0.5557 1.0000;
 17300  0.4178 0.5763 1.0000;
 17400  0.4429 0.5547 1.0000;
 17400  0.4167 0.5755 1.0000;
 17500  0.4418 0.5538 1.0000;
 17500  0.4156 0.5746 1.0000;
 17600  0.4408 0.5528 1.0000;
 17600  0.4145 0.5738 1.0000;
 17700  0.4397 0.5519 1.0000;
 17700  0.4134 0.5729 1.0000;
 17800  0.4387 0.5510 1.0000;
 17800  0.4124 0.5721 1.0000;
 17900  0.4377 0.5501 1.0000;
 17900  0.4113 0.5713 1.0000;
 18000  0.4367 0.5492 1.0000;
 18000  0.4103 0.5705 1.0000;
 18100  0.4357 0.5483 1.0000;
 18100  0.4093 0.5697 1.0000;
 18200  0.4348 0.5475 1.0000;
 18200  0.4083 0.5689 1.0000;
 18300  0.4338 0.5466 1.0000;
 18300  0.4074 0.5681 1.0000;
 18400  0.4329 0.5458 1.0000;
 18400  0.4064 0.5674 1.0000;
 18500  0.4319 0.5450 1.0000;
 18500  0.4055 0.5666 1.0000;
 18600  0.4310 0.5442 1.0000;
 18600  0.4045 0.5659 1.0000;
 18700  0.4301 0.5434 1.0000;
 18700  0.4036 0.5652 1.0000;
 18800  0.4293 0.5426 1.0000;
 18800  0.4027 0.5644 1.0000;
 18900  0.4284 0.5418 1.0000;
 18900  0.4018 0.5637 1.0000;
 19000  0.4275 0.5410 1.0000;
 19000  0.4009 0.5630 1.0000;
 19100  0.4267 0.5403 1.0000;
 19100  0.4001 0.5623 1.0000;
 19200  0.4258 0.5395 1.0000;
 19200  0.3992 0.5616 1.0000;
 19300  0.4250 0.5388 1.0000;
 19300  0.3984 0.5610 1.0000;
 19400  0.4242 0.5381 1.0000;
 19400  0.3975 0.5603 1.0000;
 19500  0.4234 0.5373 1.0000;
 19500  0.3967 0.5596 1.0000;
 19600  0.4226 0.5366 1.0000;
 19600  0.3959 0.5590 1.0000;
 19700  0.4218 0.5359 1.0000;
 19700  0.3951 0.5584 1.0000;
 19800  0.4211 0.5352 1.0000;
 19800  0.3943 0.5577 1.0000;
 19900  0.4203 0.5345 1.0000;
 19900  0.3935 0.5571 1.0000;
 20000  0.4196 0.5339 1.0000;
 20000  0.3928 0.5565 1.0000;
 20100  0.4188 0.5332 1.0000;
 20100  0.3920 0.5559 1.0000;
 20200  0.4181 0.5325 1.0000;
 20200  0.3913 0.5553 1.0000;
 20300  0.4174 0.5319 1.0000;
 20300  0.3905 0.5547 1.0000;
 20400  0.4167 0.5312 1.0000;
 20400  0.3898 0.5541 1.0000;
 20500  0.4160 0.5306 1.0000;
 20500  0.3891 0.5535 1.0000;
 20600  0.4153 0.5300 1.0000;
 20600  0.3884 0.5529 1.0000;
 20700  0.4146 0.5293 1.0000;
 20700  0.3877 0.5524 1.0000;
 20800  0.4139 0.5287 1.0000;
 20800  0.3870 0.5518 1.0000;
 20900  0.4133 0.5281 1.0000;
 20900  0.3863 0.5513 1.0000;
 21000  0.4126 0.5275 1.0000;
 21000  0.3856 0.5507 1.0000;
 21100  0.4119 0.5269 1.0000;
 21100  0.3850 0.5502 1.0000;
 21200  0.4113 0.5264 1.0000;
 21200  0.3843 0.5496 1.0000;
 21300  0.4107 0.5258 1.0000;
 21300  0.3836 0.5491 1.0000;
 21400  0.4100 0.5252 1.0000;
 21400  0.3830 0.5486 1.0000;
 21500  0.4094 0.5246 1.0000;
 21500  0.3824 0.5481 1.0000;
 21600  0.4088 0.5241 1.0000;
 21600  0.3817 0.5476 1.0000;
 21700  0.4082 0.5235 1.0000;
 21700  0.3811 0.5471 1.0000;
 21800  0.4076 0.5230 1.0000;
 21800  0.3805 0.5466 1.0000;
 21900  0.4070 0.5224 1.0000;
 21900  0.3799 0.5461 1.0000;
 22000  0.4064 0.5219 1.0000;
 22000  0.3793 0.5456 1.0000;
 22100  0.4059 0.5214 1.0000;
 22100  0.3787 0.5451 1.0000;
 22200  0.4053 0.5209 1.0000;
 22200  0.3781 0.5446 1.0000;
 22300  0.4047 0.5203 1.0000;
 22300  0.3776 0.5441 1.0000;
 22400  0.4042 0.5198 1.0000;
 22400  0.3770 0.5437 1.0000;
 22500  0.4036 0.5193 1.0000;
 22500  0.3764 0.5432 1.0000;
 22600  0.4031 0.5188 1.0000;
 22600  0.3759 0.5428 1.0000;
 22700  0.4026 0.5183 1.0000;
 22700  0.3753 0.5423 1.0000;
 22800  0.4020 0.5178 1.0000;
 22800  0.3748 0.5419 1.0000;
 22900  0.4015 0.5174 1.0000;
 22900  0.3742 0.5414 1.0000;
 23000  0.4010 0.5169 1.0000;
 23000  0.3737 0.5410 1.0000;
 23100  0.4005 0.5164 1.0000;
 23100  0.3732 0.5405 1.0000;
 23200  0.4000 0.5159 1.0000;
 23200  0.3726 0.5401 1.0000;
 23300  0.3995 0.5155 1.0000;
 23300  0.3721 0.5397 1.0000;
 23400  0.3990 0.5150 1.0000;
 23400  0.3716 0.5393 1.0000;
 23500  0.3985 0.5146 1.0000;
 23500  0.3711 0.5389 1.0000;
 23600  0.3980 0.5141 1.0000;
 23600  0.3706 0.5384 1.0000;
 23700  0.3975 0.5137 1.0000;
 23700  0.3701 0.5380 1.0000;
 23800  0.3970 0.5132 1.0000;
 23800  0.3696 0.5376 1.0000;
 23900  0.3966 0.5128 1.0000;
 23900  0.3692 0.5372 1.0000;
 24000  0.3961 0.5123 1.0000;
 24000  0.3687 0.5368 1.0000;
 24100  0.3956 0.5119 1.0000;
 24100  0.3682 0.5365 1.0000;
 24200  0.3952 0.5115 1.0000;
 24200  0.3677 0.5361 1.0000;
 24300  0.3947 0.5111 1.0000;
 24300  0.3673 0.5357 1.0000;
 24400  0.3943 0.5107 1.0000;
 24400  0.3668 0.5353 1.0000;
 24500  0.3938 0.5103 1.0000;
 24500  0.3664 0.5349 1.0000;
 24600  0.3934 0.5098 1.0000;
 24600  0.3659 0.5346 1.0000;
 24700  0.3930 0.5094 1.0000;
 24700  0.3655 0.5342 1.0000;
 24800  0.3925 0.5090 1.0000;
 24800  0.3650 0.5338 1.0000;
 24900  0.3921 0.5086 1.0000;
 24900  0.3646 0.5335 1.0000;
 25000  0.3917 0.5083 1.0000;
 25000  0.3642 0.5331 1.0000;
 25100  0.3913 0.5079 1.0000;
 25100  0.3637 0.5328 1.0000;
 25200  0.3909 0.5075 1.0000;
 25200  0.3633 0.5324 1.0000;
 25300  0.3905 0.5071 1.0000;
 25300  0.3629 0.5321 1.0000;
 25400  0.3901 0.5067 1.0000;
 25400  0.3625 0.5317 1.0000;
 25500  0.3897 0.5064 1.0000;
 25500  0.3621 0.5314 1.0000;
 25600  0.3893 0.5060 1.0000;
 25600  0.3617 0.5310 1.0000;
 25700  0.3889 0.5056 1.0000;
 25700  0.3613 0.5307 1.0000;
 25800  0.3885 0.5053 1.0000;
 25800  0.3609 0.5304 1.0000;
 25900  0.3881 0.5049 1.0000;
 25900  0.3605 0.5300 1.0000;
 26000  0.3877 0.5045 1.0000;
 26000  0.3601 0.5297 1.0000;
 26100  0.3874 0.5042 1.0000;
 26100  0.3597 0.5294 1.0000;
 26200  0.3870 0.5038 1.0000;
 26200  0.3593 0.5291 1.0000;
 26300  0.3866 0.5035 1.0000;
 26300  0.3589 0.5288 1.0000;
 26400  0.3863 0.5032 1.0000;
 26400  0.3586 0.5284 1.0000;
 26500  0.3859 0.5028 1.0000;
 26500  0.3582 0.5281 1.0000;
 26600  0.3855 0.5025 1.0000;
 26600  0.3578 0.5278 1.0000;
 26700  0.3852 0.5021 1.0000;
 26700  0.3575 0.5275 1.0000;
 26800  0.3848 0.5018 1.0000;
 26800  0.3571 0.5272 1.0000;
 26900  0.3845 0.5015 1.0000;
 26900  0.3567 0.5269 1.0000;

 27000  0.3841 0.5012 1.0000;
 27000  0.3564 0.5266 1.0000;
 27100  0.3838 0.5008 1.0000;
 27100  0.3560 0.5263 1.0000;
 27200  0.3835 0.5005 1.0000;
 27200  0.3557 0.5260 1.0000;
 27300  0.3831 0.5002 1.0000;
 27300  0.3553 0.5257 1.0000;
 27400  0.3828 0.4999 1.0000;
 27400  0.3550 0.5255 1.0000;
 27500  0.3825 0.4996 1.0000;
 27500  0.3546 0.5252 1.0000;
 27600  0.3821 0.4993 1.0000;
 27600  0.3543 0.5249 1.0000;
 27700  0.3818 0.4990 1.0000;
 27700  0.3540 0.5246 1.0000;
 27800  0.3815 0.4987 1.0000;
 27800  0.3536 0.5243 1.0000;
 27900  0.3812 0.4984 1.0000;
 27900  0.3533 0.5241 1.0000;
 28000  0.3809 0.4981 1.0000;
 28000  0.3530 0.5238 1.0000;
 28100  0.3805 0.4978 1.0000;
 28100  0.3527 0.5235 1.0000;
 28200  0.3802 0.4975 1.0000;
 28200  0.3524 0.5232 1.0000;
 28300  0.3799 0.4972 1.0000;
 28300  0.3520 0.5230 1.0000;
 28400  0.3796 0.4969 1.0000;
 28400  0.3517 0.5227 1.0000;
 28500  0.3793 0.4966 1.0000;
 28500  0.3514 0.5225 1.0000;
 28600  0.3790 0.4963 1.0000;
 28600  0.3511 0.5222 1.0000;
 28700  0.3787 0.4960 1.0000;
 28700  0.3508 0.5219 1.0000;
 28800  0.3784 0.4958 1.0000;
 28800  0.3505 0.5217 1.0000;
 28900  0.3781 0.4955 1.0000;
 28900  0.3502 0.5214 1.0000;
 29000  0.3779 0.4952 1.0000;
 29000  0.3499 0.5212 1.0000;
 29100  0.3776 0.4949 1.0000;
 29100  0.3496 0.5209 1.0000;
 29200  0.3773 0.4947 1.0000;
 29200  0.3493 0.5207 1.0000;
 29300  0.3770 0.4944 1.0000;
 29300  0.3490 0.5204 1.0000;
 29400  0.3767 0.4941 1.0000;
 29400  0.3487 0.5202 1.0000;
 29500  0.3764 0.4939 1.0000;
 29500  0.3485 0.5200 1.0000;
 29600  0.3762 0.4936 1.0000;
 29600  0.3482 0.5197 1.0000;
 29700  0.3759 0.4934 1.0000;
 29700  0.3479 0.5195 1.0000;
 29800  0.3756 0.4931 1.0000;
 29800  0.3476 0.5192 1.0000;
 29900  0.3754 0.4928 1.0000;
 29900  0.3473 0.5190 1.0000;
 30000  0.3751 0.4926 1.0000;
 30000  0.3471 0.5188 1.0000;
 30100  0.3748 0.4923 1.0000;
 30100  0.3468 0.5186 1.0000;
 30200  0.3746 0.4921 1.0000;
 30200  0.3465 0.5183 1.0000;
 30300  0.3743 0.4918 1.0000;
 30300  0.3463 0.5181 1.0000;
 30400  0.3741 0.4916 1.0000;
 30400  0.3460 0.5179 1.0000;
 30500  0.3738 0.4914 1.0000;
 30500  0.3457 0.5177 1.0000;
 30600  0.3735 0.4911 1.0000;
 30600  0.3455 0.5174 1.0000;
 30700  0.3733 0.4909 1.0000;
 30700  0.3452 0.5172 1.0000;
 30800  0.3730 0.4906 1.0000;
 30800  0.3450 0.5170 1.0000;
 30900  0.3728 0.4904 1.0000;
 30900  0.3447 0.5168 1.0000;
 31000  0.3726 0.4902 1.0000;
 31000  0.3444 0.5166 1.0000;
 31100  0.3723 0.4899 1.0000;
 31100  0.3442 0.5164 1.0000;
 31200  0.3721 0.4897 1.0000;
 31200  0.3439 0.5161 1.0000;
 31300  0.3718 0.4895 1.0000;
 31300  0.3437 0.5159 1.0000;
 31400  0.3716 0.4893 1.0000;
 31400  0.3435 0.5157 1.0000;
 31500  0.3714 0.4890 1.0000;
 31500  0.3432 0.5155 1.0000;
 31600  0.3711 0.4888 1.0000;
 31600  0.3430 0.5153 1.0000;
 31700  0.3709 0.4886 1.0000;
 31700  0.3427 0.5151 1.0000;
 31800  0.3707 0.4884 1.0000;
 31800  0.3425 0.5149 1.0000;
 31900  0.3704 0.4881 1.0000;
 31900  0.3423 0.5147 1.0000;
 32000  0.3702 0.4879 1.0000;
 32000  0.3420 0.5145 1.0000;
 32100  0.3700 0.4877 1.0000;
 32100  0.3418 0.5143 1.0000;
 32200  0.3698 0.4875 1.0000;
 32200  0.3416 0.5141 1.0000;
 32300  0.3695 0.4873 1.0000;
 32300  0.3413 0.5139 1.0000;
 32400  0.3693 0.4871 1.0000;
 32400  0.3411 0.5137 1.0000;
 32500  0.3691 0.4869 1.0000;
 32500  0.3409 0.5135 1.0000;
 32600  0.3689 0.4867 1.0000;
 32600  0.3407 0.5133 1.0000;
 32700  0.3687 0.4864 1.0000;
 32700  0.3404 0.5132 1.0000;
 32800  0.3684 0.4862 1.0000;
 32800  0.3402 0.5130 1.0000;
 32900  0.3682 0.4860 1.0000;
 32900  0.3400 0.5128 1.0000;
 33000  0.3680 0.4858 1.0000;
 33000  0.3398 0.5126 1.0000;
 33100  0.3678 0.4856 1.0000;
 33100  0.3396 0.5124 1.0000;
 33200  0.3676 0.4854 1.0000;
 33200  0.3393 0.5122 1.0000;
 33300  0.3674 0.4852 1.0000;
 33300  0.3391 0.5120 1.0000;
 33400  0.3672 0.4850 1.0000;
 33400  0.3389 0.5119 1.0000;
 33500  0.3670 0.4848 1.0000;
 33500  0.3387 0.5117 1.0000;
 33600  0.3668 0.4847 1.0000;
 33600  0.3385 0.5115 1.0000;
 33700  0.3666 0.4845 1.0000;
 33700  0.3383 0.5113 1.0000;
 33800  0.3664 0.4843 1.0000;
 33800  0.3381 0.5112 1.0000;
 33900  0.3662 0.4841 1.0000;
 33900  0.3379 0.5110 1.0000;
 34000  0.3660 0.4839 1.0000;
 34000  0.3377 0.5108 1.0000;
 34100  0.3658 0.4837 1.0000;
 34100  0.3375 0.5106 1.0000;
 34200  0.3656 0.4835 1.0000;
 34200  0.3373 0.5105 1.0000;
 34300  0.3654 0.4833 1.0000;
 34300  0.3371 0.5103 1.0000;
 34400  0.3652 0.4831 1.0000;
 34400  0.3369 0.5101 1.0000;
 34500  0.3650 0.4830 1.0000;
 34500  0.3367 0.5100 1.0000;
 34600  0.3649 0.4828 1.0000;
 34600  0.3365 0.5098 1.0000;
 34700  0.3647 0.4826 1.0000;
 34700  0.3363 0.5096 1.0000;
 34800  0.3645 0.4824 1.0000;
 34800  0.3361 0.5095 1.0000;
 34900  0.3643 0.4822 1.0000;
 34900  0.3359 0.5093 1.0000;
 35000  0.3641 0.4821 1.0000;
 35000  0.3357 0.5091 1.0000;
 35100  0.3639 0.4819 1.0000;
 35100  0.3356 0.5090 1.0000;
 35200  0.3638 0.4817 1.0000;
 35200  0.3354 0.5088 1.0000;
 35300  0.3636 0.4815 1.0000;
 35300  0.3352 0.5087 1.0000;
 35400  0.3634 0.4814 1.0000;
 35400  0.3350 0.5085 1.0000;
 35500  0.3632 0.4812 1.0000;
 35500  0.3348 0.5084 1.0000;
 35600  0.3630 0.4810 1.0000;
 35600  0.3346 0.5082 1.0000;
 35700  0.3629 0.4809 1.0000;
 35700  0.3345 0.5080 1.0000;
 35800  0.3627 0.4807 1.0000;
 35800  0.3343 0.5079 1.0000;
 35900  0.3625 0.4805 1.0000;
 35900  0.3341 0.5077 1.0000;
 36000  0.3624 0.4804 1.0000;
 36000  0.3339 0.5076 1.0000;
 36100  0.3622 0.4802 1.0000;
 36100  0.3338 0.5074 1.0000;
 36200  0.3620 0.4800 1.0000;
 36200  0.3336 0.5073 1.0000;
 36300  0.3619 0.4799 1.0000;
 36300  0.3334 0.5071 1.0000;
 36400  0.3617 0.4797 1.0000;
 36400  0.3332 0.5070 1.0000;
 36500  0.3615 0.4796 1.0000;
 36500  0.3331 0.5068 1.0000;
 36600  0.3614 0.4794 1.0000;
 36600  0.3329 0.5067 1.0000;
 36700  0.3612 0.4792 1.0000;
 36700  0.3327 0.5066 1.0000;
 36800  0.3610 0.4791 1.0000;
 36800  0.3326 0.5064 1.0000;
 36900  0.3609 0.4789 1.0000;
 36900  0.3324 0.5063 1.0000;
 37000  0.3607 0.4788 1.0000;
 37000  0.3322 0.5061 1.0000;
 37100  0.3605 0.4786 1.0000;
 37100  0.3321 0.5060 1.0000;
 37200  0.3604 0.4785 1.0000;
 37200  0.3319 0.5058 1.0000;
 37300  0.3602 0.4783 1.0000;
 37300  0.3317 0.5057 1.0000;
 37400  0.3601 0.4782 1.0000;
 37400  0.3316 0.5056 1.0000;
 37500  0.3599 0.4780 1.0000;
 37500  0.3314 0.5054 1.0000;
 37600  0.3598 0.4779 1.0000;
 37600  0.3313 0.5053 1.0000;
 37700  0.3596 0.4777 1.0000;
 37700  0.3311 0.5052 1.0000;
 37800  0.3595 0.4776 1.0000;
 37800  0.3309 0.5050 1.0000;
 37900  0.3593 0.4774 1.0000;
 37900  0.3308 0.5049 1.0000;
 38000  0.3592 0.4773 1.0000;
 38000  0.3306 0.5048 1.0000;
 38100  0.3590 0.4771 1.0000;
 38100  0.3305 0.5046 1.0000;
 38200  0.3589 0.4770 1.0000;
 38200  0.3303 0.5045 1.0000;
 38300  0.3587 0.4768 1.0000;
 38300  0.3302 0.5044 1.0000;
 38400  0.3586 0.4767 1.0000;
 38400  0.3300 0.5042 1.0000;
 38500  0.3584 0.4766 1.0000;
 38500  0.3299 0.5041 1.0000;
 38600  0.3583 0.4764 1.0000;
 38600  0.3297 0.5040 1.0000;
 38700  0.3581 0.4763 1.0000;
 38700  0.3296 0.5038 1.0000;
 38800  0.3580 0.4761 1.0000;
 38800  0.3294 0.5037 1.0000;
 38900  0.3579 0.4760 1.0000;
 38900  0.3293 0.5036 1.0000;
 39000  0.3577 0.4759 1.0000;
 39000  0.3291 0.5035 1.0000;
 39100  0.3576 0.4757 1.0000;
 39100  0.3290 0.5033 1.0000;
 39200  0.3574 0.4756 1.0000;
 39200  0.3288 0.5032 1.0000;
 39300  0.3573 0.4755 1.0000;
 39300  0.3287 0.5031 1.0000;
 39400  0.3572 0.4753 1.0000;
 39400  0.3286 0.5030 1.0000;
 39500  0.3570 0.4752 1.0000;
 39500  0.3284 0.5028 1.0000;
 39600  0.3569 0.4751 1.0000;
 39600  0.3283 0.5027 1.0000;
 39700  0.3567 0.4749 1.0000;
 39700  0.3281 0.5026 1.0000;
 39800  0.3566 0.4748 1.0000;
 39800  0.3280 0.5025 1.0000;
 39900  0.3565 0.4747 1.0000;
 39900  0.3279 0.5024 1.0000;
 40000  0.3563 0.4745 1.0000;
 40000  0.3277 0.5022 1.0000];
