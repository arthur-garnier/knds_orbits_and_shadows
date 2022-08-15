User manual for the package "knds_orbits_and_shadows".

This is a package for Scilab 6.1.1, loaded with IPCV 4.1.2, to draw the shadow of a Kerr--Newman--(anti) de Sitter (KNdS) black hole, possibly equipped with a thin Keplerian accretion disk, radiating as a blackbody. The shadow can be either drawn from any standard image, or only the accretion disk is drawn on a black background.
This package also allows to draw (massive or null) orbits in a KNdS space-time, using different integration methods of the geodesic equation.

---------------------------------------------------------------------------------------------------

First of all, install Scilab 6.1.1 and the package IPCV 4.1.2 (https://ipcv.scilab-academy.com/).
On Windows, just install Scilab 6.1.1 via https://www.scilab.org/download/scilab-6.1.1 and install IPCV using the ATOMS Module Manager. Re-launch Scilab and everything should be OK.
This package can be a bit tough to install on Linux x64, but once IPCV has been installed (using ATOMS), the debugging instructions and patches by the author can be found at https://atoms.scilab.org/toolboxes/IPCV.
When launching Scilab on a terminal, once in the Scilab installation directory, enter the instruction "LIBGL_ALWAYS_SOFTWARE=1 ./bin/scilab".
This ensures that the graphics done using IPCV will indeed be displayed.

Next, put the content of the present folder anywhere and in the example.sce file, change the first line to match the directory of the files (and images!).
Execute the file example.sce; it displays all the functions of the programs, so it should be a good indicator of the sanity of the package. The full execution takes about 35 seconds on a 12-core 2.60 GHz CPU with 16 Go of RAM.

---------------------------------------------------------------------------------------------------

Description of the main functions of the package:


- First, the file 'aux.sci' contains all the auxiliary functions that are needed for the computations, namely the metric matrices, their derivatives and the Christoffel symbols.
It also contains the library for the black-body radiation colors from http://www.vendian.org/mncharity/dir3/blackbody/.


- orbit.sci computes the trajectory of a test particle in a KNdS space-time.

The synthax is as follows: [Vecc,HAM,CAR]=orbit(Lambda,Mass,Kerr,Newman,IniConds,Form,Tau,N,Mu,Conserv), where Lambda is the cosmological constant, Mass is the mass of the black hole, Kerr is its Kerr parameter and Newman its charge.
The vector IniConds records the initial datum of the geodesic, in Boyer-Lindquist coordinates and SI units, IniConds=(r0,theta0,phi0,\dot{r}0,\dot{theta}0,\dot{\phi}0).
The variable Form denotes the formulation to take for the computation (a string with value "Polar", "Weierstrass" (for Lambda=0), "Euler-Lagrange", "Carter", "Hamilton", "Symplectic Euler p", "Symplectic Euler q", "Verlet" or "Stormer-Verlet").
The variable Tau denotes the maximal affine parameter where the trajectory is computed, N is the number of discretization points, Mu is the "massiveness" of the particle: Mu=1 for a massive particle and Mu=0 for a photon.
Finally, the variable Conserv is set to 1 if the user is willing to compute the Hamiltonian and Carter constant at each node of the discretization. Otherwise, the user is invited to set Conserv to 0.

The output is as follows: the vector Vecc contains the position (r,theta,phi) (in SI units) of the trajectory at each node k*Tau/N (0<k<N) of the discretization. 
If Conserv is set to 1, then the vectors HAM and CAR respectively contain the value of the Hamiltonian and Carter constant at each node.


- shadow.sci shadows a KNdS black hole, with a standard image (jpeg, png...), or an accretion disk, or both.

The synthax is shadow(Lambda,Mass,Kerr,Newman,Image,Accretion_data), where Lambda is the cosmological constant, M is the mass, Kerr is the Kerr parameter and Newman is the charge.
The variable Image is a string formed with the name (with extension) of the picture to transform (the file should be in the same folder as the functions).
The variable Accretion_data is a list with seven entries.
Its first entry is a non-negative integer: set it to 0 yields the shadow without any accretion disk, set it to 1 gives the picture with the accretion disk and otherwise, the value should be even and only the shadow of the accretion disk is drawn, with a resolution equal to the chosen even integer.
The second entry of the list is the inclination angle from the equatorial plane (so that set this angle to 0 means to shadow the black hole, as seen from the equatorial plane).
The third entry is a string setting the type of radiation required. Set it to " " only computes the effects (graviational, Doppler, both, see below) and yields the shift along the disk, displayed as a shade of colors from red to blue.
    Set this variable to "Black-body" computes the temperature as a blackbody radiation. Finally, set it to "Custom" allows to specify the inner and outer temperatures (see below).
The fourth variable is a string too, specifying the various shifts to take into account. Set it to "Gravitation" (resp. to "Doppler") only computes the gravitational (resp. Doppler) shift. To take both effects into account, set this variable to "Doppler+". If " " is chosen, none of these effects is computed.
    In any other case than the ones just described for the third and/or fourth variable, the color of the accretion disk is arbitrarily set to [R,G,B]=[255,69,0] and the brightness is computed with a linear scale from the outer to the inner radius.
The fifth variable is a vector with two entries: the respective inner and outer radii of the disk (in terms of the Schwarzschild radius).
The sixth variable is a vector with one or three entries: the first one is the accretion rate and the other ones are (if specified) the respective inner and outer temperature of the disk (in Kelvin). These two values are needed only in the case where the option "Custom" is chosen and are ignored otherwise.
Finally, the seventh variable is a non-negative integer: the brightness scaling. If it is set to 0, then the brightness is linearly computed as above. Otherwise, it is computed with Planck's law and rescaled using this value. This is to be adjusted case by case.

The output is just the computed picture, in a Scilab figure.


- shadow_full.sci and shadow_wp.sci

The synthax, input and output are exactly the same as for shadow.sci.
The difference is that the function shadow.sci doesn't use the Weierstrass function to simplify the computation of the coordinate range of the pixels. So this program is more efficient for images with a very low resolution, but it is much slower than the function shadow.sci.
The function shadow_wp.sci only uses the Weierstrass function and is thus well suited for non-rotating (a=0) black holes only. The user is encouraged to prefer this function for non-rotating black holes, as it is much faster than shadow.sci.

---------------------------------------------------------------------------------------------------

For more details on the equations and modelization, the reader is refered to the article available at {ARXIV URL}.
For any question, suggestion, commentary, remark, the user is invited to contact the author by email at art-uhr[at]orange[dot]fr.
