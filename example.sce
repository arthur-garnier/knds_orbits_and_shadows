//Before all, change the current directory to the one containing the scripts and execute the functions:
chdir("/home/arthur/Documents/BH/GitHub");//The absolute path to the package directory should be put here: chdir("PATH");
exec('auxi.sci', -1); exec('orbit.sci', -1); exec('shadow.sci', -1);

//We test the functions for Lambda=0 and Lambda<>0
for Lambda=[0,3.3e-4]
    //First, we define the required parameters and initial data:
    cSI=299792458; GSI=6.67408e-11; M=4e30; Rs=2*GSI*M/cSI^2; a=0.95; Q=0.3;

    //The initial value of the (massive) orbit:
    X=[30656,%pi/2,0,0,1677.2693,2234.125]'; tau=0.0075; N=1000; mu=1;

    //Different formulations (with colors to distinguish between them on the graphs):
    eqnss=["Euler-Lagrange","Hamilton","Carter","Verlet","Stormer-Verlet","Symplectic Euler p","Symplectic Euler q"];
    colors=["blue","green","red","cyan","purple","orange","black","magenta","pink"];

    //The Schwarzszchild sphere:
    deff("[x,y,z]=sph(alp,tet)",["x=Rs*cos(alp).*cos(tet)";"y=Rs*cos(alp).*sin(tet)";"z=Rs*sin(alp)"]); [xx,yy,zz]=eval3dp(sph,linspace(-%pi/2,%pi/2,40),linspace(0,%pi*2,20));

    //On a new figure, draw the Schwarzschild sphere and the orbit using the various formulations:
    figure(); plot3d(xx,yy,zz); times=[]; HAMS=[]; CARS=[];
    for i=[1:size(eqnss)(2)]
        tic(); [Vec,HAM,CAR]=orbit(Lambda,M,a,Q,X,eqnss(i),tau,N,mu,1,0); ttt=toc(); times=[times,ttt];
        R=Vec(1,:); theta=Vec(2,:); phi=Vec(3,:); Mm=[1;X(4);X(5);X(6)]; HAMS(:,$+1)=HAM/HAM(2); CARS(:,$+1)=CAR/CAR(2);
        param3d(sqrt(R.^2+a^2).*sin(theta).*cos(phi),sqrt(R.^2+a^2).*sin(theta).*sin(phi),R.*cos(theta),"x@y@z"); curve = gce(); curve.foreground = color(colors(i));
    end
    legend(eqnss,1); xtitle('Comparison between different orbits (mu='+mtlb_num2str(-mu)+')','X','Y','Z');
    aa=gca(); aa.isoview="on"; aa = gcf(); aa.background = 8; aa=gce(); aa=aa.children(1); aa.background=8; aa.font_color=-1;

    //Compare the conservation of the Hamiltonian:
    figure(); I=cSI*2/Rs*[tau/N:tau/N:tau]';
    plot2d(I,HAMS(2:$,:),style=[color(colors(1)),color(colors(2)),color(colors(3)),color(colors(4)),color(colors(5)),color(colors(6)),color(colors(7)),color(colors(8)),color(colors(9))]);
    legend(eqnss,1);
    xtitle('Hamiltonian conservation for different methods (mu='+mtlb_num2str(-mu)+')','proper time[M/c]','H/H0 where H=1/2 g^ij p_i p_j');
    aa = gcf(); aa.background = 8; aa=gce(); aa=aa.children(1); aa.background=8; aa.font_color=-1;

    //Compare the conservation of the Carter constant:
    figure(); I=cSI*2/Rs*[tau/N:tau/N:tau]';
    plot2d(I,CARS(2:$,:),style=[color(colors(1)),color(colors(2)),color(colors(3)),color(colors(4)),color(colors(5)),color(colors(6)),color(colors(7)),color(colors(8)),color(colors(9))]);
    legend(eqnss,1);
    xtitle('Carter constant conservation for different methods (mu='+mtlb_num2str(-mu)+')','proper time[M/c]','C/C0 where C is the Carter constant');
    aa = gcf(); aa.background = 8; aa=gce(); aa=aa.children(1); aa.background=8; aa.font_color=-1;

    //Maximal deviation (Hamilton and Carter) for each method:
    dev_ham=[]; dev_car=[];
    for i=[1:size(eqnss)(2)]
        dev_ham=[dev_ham,max(abs(HAMS(2:$,i)-1))]; dev_car=[dev_car,max(abs(CARS(2:$,i)-1))];
    end

    //------------------------------------------------------------------------------

    //Testing the shadowing programs (the reader is invited to un-comment the three lines below to test the effect of the two shifts we introduced):
    Mass=4e30; Kerr=0.95; Newman=0.3; Image='figure24.png';
    Accretion=list(1,%pi/18,"Black-body","Doppler+",[1.455,6],[3],3800,1);
    //Accretion=list(16,%pi/18," ","Gravitation",[1.455,6],[3],3800,0);
    //Accretion=list(16,%pi/18," ","Doppler",[1.455,6],[3],3800,0);
    //Accretion=list(16,%pi/18," ","Doppler+",[1.455,6],[3],3800,0);
    tic(); shadow(Lambda,Mass,Kerr,Newman,Image,Accretion); t_full=toc();
    //tic(); shadow(Lambda,Mass,Kerr,Newman,Image,Accretion); t=toc();
    tic(); shadow(Lambda,Mass,0,Newman,Image,Accretion); t_wp=toc();
end
