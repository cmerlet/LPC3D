/* This program is a lattice model to calculate quantities of adsorbed ions, diffusion coefficients and predict NMR spectra of ions adsorbed in 
porous carbon matrices:
 - for diffusion coefficients calculations, data derived from molecular dynamics simulations is used as input to determine in-pore adsorption profiles, 
   and data from experiments is used to obtain energy barriers governing transitions between lattice sites;
 - for NMR spectra calculations, the model uses chemical shifts obtained from density functional theory calculations.
In this model, only one carbon electrode is simulated and its structure is represented as a cubic tridimensional set of inter-connected
discrete sites, separated by a lattice spacing, "a". */

/* Information describing the lattice model for porous carbons is related here:
J. Chem. Phys., 142, 094701 (2015)
and available at:
http://arxiv.org/abs/1412.7043 
and 
http://scitation.aip.org/content/aip/journal/jcp/142/9/10.1063/1.4913368 */

/* A description of the model with applied potentials is related here:
Electrochimica Acta, 327, 135022 (2019)
and available at:
https://arxiv.org/abs/1910.02663
and
https://www.sciencedirect.com/science/article/pii/S0013468619318936 */

/* Since March 2023, it is possible to include a bulk region around the carbon particle. */

/* To execute the program, you have to launch it with an input file (./program < inputfile).
The beginning of the file is the same for all types of simulations:
1st line: numbers of lattice sites in x, y and z directions (> Nx, Ny, Nz)
2nd line: lattice parameter in anstroms (> a)
3nd line: dwell time in seconds (> dwell)
4th line: number of steps between measures (> nsample)
5th line: temperature in Kelvin (> T)
6th line: Larmor frequency of the studied nucleus under the magnetic field considered in MHz (> larmorfreq)
This frequency can be found on http://nmr.cemhti.cnrs-orleans.fr/dmfit/Tables/TableIsotopes.aspx.
7th line: random seed for the steps involving random (> seed)
8th line: obstacles (> test_latt)
0: no obstacles, all sites are accessible
1: lattice given from a previous run or built differently
2: obstacles assigned randomly, some sites are unaccessible
3: obstacles defined according to positions of carbon atoms
9th line: no 9th line if test_latt=0, name of the obstacles' file if test_latt=1 (> namelatt), 
coverage if test_latt=2 (> coverage), name of the positions' file and minimum distance if 
test_latt=3 (> namepos, Dmin)
(positions' file = xyz file with coordinates in the same units as the lattice parameter)
10th line: number of values for the Fourier transform (> Nvalues), please use an even number of values
11th line: pulse sequence (> pseq)
12th line: type of the simulation, see below (> simtype)
*/

/* At the moment, the pulse sequence can be:
- "single": single pulse;
- "echo": spin echo = 90º - tau - 180º - tau 
AT THIS DATE THE ECHO PULSE SEQUENCE HAS NOT BEEN VALIDATED! */

/* FROM THIS POINT, THE FILES DIFFER DEPENDING ON THE TYPE OF SIMULATION */
/* The types are the following:
- "planar": one site = one position in a pore, note that this simulation type assigns chemical shift depending on 
the z coordinate only;
- "porous": one site = one position in a pore, at this point, the chemical shift is assigned depending on the 
smallest 'considered site' - 'non accessible site' distance;
- "coarse": one site = one pore. */


/************************* Planar *************************/
/* 
13th line: name of the file with the free energies (> nameeng), 1st column should be in anstroms, 2nd column should be in kJ/mol 
14th line: number of free energies values (> Neng)
15th line: do we use the values directly, do we build a pore or do we use Xing fit (Carbon, 77, 1132, 2014)? (> poreshift) 
16th line: name of the file with the frequencies if direct/indirect before (> namefreq1), 1st column should be in anstroms
17th line: number of frequencies if direct/indirect before (> Nfreq)
*/ 

/************************* Porous *************************/
/*
13th line: type of free energies assignment (> test_feng)
 test_feng can be "1D" if we use a one dimensional profile or "3D" if we use a three dimensional file
 coordinates should in anstroms, free energies should be in kJ/mol 
14th line: name of the free energies file and number of values (> nameng, Neng)
15th line: type of chemical shifts (> test_shift)
 test_shift can be "1D" if we use a one dimensional profile (for example from Gaussian), "Xing" if we use Xing model, 
 "3D" if we use a three dimensional file, "dipolar" if we use the dipolar model
16th line: name of the file with the chemical shifts and number of values if "1D" or "3D" before (> namefreq1, Nfreq)
 coordinates should be in anstroms, chemical shifts should be in ppm
*/ 

/************************* Coarse *************************/
/*
13th line: number of integrated densities values (> Npores), it's equal to the number of frequencies
following line: how to read the integrated densities (> densread) (= distinguished, to read anion, cation and solvent separately), (= undistinguished, to read density as a whole)
following line (3): name of the file(s) with the integrated densities of anions and cations (if distinguished reading) or considered species (if undistinguished reading), unit does not really matter as long as ratios of densities are correct
(if distinguished reading) following line: Solvent consideration "solvent_ON" or "solvent_OFF" (> solvent_activ)
(if distinguished reading) following line: name of the file with the integrated densities of solvent if solvent_ON
(if distinguished reading) following line: type of molecule (> iontype = anion, cation or solvent) for which diffusion and NMR calculations will be performed 
following lines (3): names of the files with the frequencies and area of the corresponding molecules in anstroms (> namefreq1 and area[1], namefreq2 and area[2], namefreq3 and area[3])
ATTENTION THEY SHOULD BE GIVEN IN THE ORDER OF INCREASING AREAS
following line: log normal distribution or discrete values for the pore size distribution (> psdtype)
(if psdtype = "lognormal") following line: mean and standard deviation for the log normal distribution of pore sizes (> meanpsd, stdpsd)  
(if psdtype = "discrete") following line: psd reading (> psdread = "mono" for one psd (carbon type) or "mixed" for two psds (carbon types))
(if psdread = "mono") following line: name of the file with the pore size distribution and number of lines (> namepsd, nlinespsd)
(if psdread = "mixed") following line: name of the the first file with the pore size distribution and number of lines (> namepsd, nlinespsd) 
(if psdread = "mixed") following line: name of the the second file with the pore size distribution and number of lines (> namepsdm, nlinespsdm)
(if psdread = "mixed") following line: mixing ratio (> mixratio) which should give the percentage of carbon type 1 in the total system (carbon type 1 corresponds to the first declared psd)
following line: log normal distribution or discrete values for the surface distribution (> surftype)
following line: mean and standard deviation for the log normal distribution of surface sizes (> meansurf, stdsurf) if surftype = "lognormal" 
or name of the file with the surface size distribution and number of lines (> namesurf, nlinessurf) if surftype = "discrete" before
following line: type of energy barriers (> Enbar_mode) 0 ---> Random gaussian distribution, 1 ---> Energies calculated from fit of experimental data (see Forse el al., Nature Ener., 2, 16216, 2016)
following line:	parameters for the energy barriers
(if Enbar_mode = 0) mean and standard deviation for the Gaussian distribution of barrier heights (> meanEa, stdEa), 
(if Enbar_mode = 1) reference temperature at which the energy barriers are calculated (> T_ref) 
following line: gradient of pore sizes? 0 if no, 1 if yes (> gradient) 
N.B.: SIMULATIONS OF LATTICES WITH GRADIENT PORE DISTRIBUTIONS MAY PRODUCE SOME INSTABILITIES!!!
(if distinguished reading) following line: Volumes of the anion, cation and solvent molecules in this order (> Vanion, Vcation and Vsolvent) and cell voltage (> volt) if solvent is considered (solvent_ON) otherwise reading only Vanion, Vcation and volt
(if undistinguished reading) following line: Volume of the considered particule (> Vspec), and cell voltage (> volt)
(if distinguished reading) following line: molecular mass of the anion, cation and solvent in this order (> anionmass, cationmass, solventmass) in u unit (for solvent consideration: same as the reading of volumes)
(if undistinguished reading) following line: molecular mass of the particule (> mspec)
following line: number of steps for the equilibration process (> step) and equilibration threshold (> thrsld)
following line: writing density files during equilibration (> wrtdns): = wrtdns_on --> write files, = wrtdns_off --> do not write files, and number of steps between two density checks (> nsmple): = number of steps if wrtdns_on (non-zero value), = 0 if wrtdns_off
following line: minimum in-pore density to be considered for calculations (> rhotrshld). 
N.B.: To compare results of different calculations, the same minimum in-pore density should be considered for all these calculations as the chosen threshold directly impacts the considered PSD (i.e. the carbon structure) 
following line: distribution preferences of obstacles (> Cover_pref): = small --> putting obstacles on small pores, = large --> putting obstacles on large pores, = none --> putting obstacles randomly with no consideration of pore dimensions. 
N.B.: THIS OPTION SHOULD BE SPECIFIED ONLY IF test_latt=2 
following line: particle size in lattice units (> particlesize): -1 or size > lattice size --> the particle is the full box, 0 --> only bulk, between 0 and Nx --> particle and bulk 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define L 1000
#define kB 8.31036e-3 
#define Nav 6.022140857e23   /* in mol-1 */
#define graphene_surf 1336.0e20 /* in ang2/g */   
/* kB is in fact kB*0.001*Na to have it in kJ.mol-1.K-1 */ 
#define PI 3.1415926536
/* Molar magnetic susceptibility of circumcoronene from 
N. Ganguli and K. S. Krichnan, Proc. R. Soc. Lond. A, 177, 168-182 (1941),
in m3/mol */
/* TO REPLACE WITH SOMETHING CLOSE TO CIRCUMCORONENE */
#define Chipara -3.2e-9
#define Chiperp -7.5e-11
#define nbins 1000
#define Rcut 1.85
#define nzones 72
#define cutoffshift 20.0
#define Lmax 9

FILE *outspinecho;
char namelatt[50],namedens[50],nameeng[50],namefreq1[100],namefreq2[100],namefreq3[100],poreshift[50];
char namepsd[50],namesurf[50],nameecho[50],test_feng[50],test_shift[50];
char namepos[50],simtype[10],pseq[10],psdtype[10],surftype[10],wrtdns[10];
int ***Lattice,Nx,Ny,Nz,dim,dir,test_latt,seed;
int **Neigh,**Neigh_modif,*current_path,Nrings,*Histo_rings,**PATHS,*Histo_dist;
int Nl,Nf,Nmin,Npores,nsample,Nvalues,Neng,Nfreq,tau,natoms,nsmple;
double a,dwell,dtime,T,larmorfreq,coverage,sumvacf,sumvacfx,sumvacfy,sumvacfz,Dmin,Dzero;
double lattvacf,lattvacfx,lattvacfy,lattvacfz,v0,dmax,Lx,Ly,Lz;
double meanpsd,stdpsd,meansurf,stdsurf,meanEa,stdEa,int_tau0,int_tau,particlesize,fixshift;
double *area,*X,*Y,*Z,*xcom,*ycom,*zcom,*xnorm,*ynorm,*znorm;
double *xvacf,*xvacfx,*xvacfy,*xvacfz,*xsumvacf,*xsumvacfx,*xsumvacfy,*xsumvacfz;
double *yvacf,*yvacfx,*yvacfy,*yvacfz,*ysumvacf,*ysumvacfx,*ysumvacfy,*ysumvacfz;
double *zvacf,*zvacfx,*zvacfy,*zvacfz,*zsumvacf,*zsumvacfx,*zsumvacfy,*zsumvacfz;
double ***Esite,**intdens,**freq,**free_eng;
double ***wij,****ImG,****ReG,*ImGtot,*ReGtot,*ImGtotsmp,*ReGtotsmp;
double *ReGcf,*ImGcf,*ReGsmpcf,*ImGsmpcf;
double ****dens,*****probv0,****alpha;

int Enbar_mode,step;
double Vanion,Vcation,Vsolvent;
double T_ref;
double ***sitesize,***pvol,***vol_occ,***surface_pore;
double ***Ntot,***pVtot,***inpore_c,***pinpore_c;
double *msdvacf,*msdvacfx,*msdvacfy,*msdvacfz;
double ****rho;
double ****dens1,****dens2;
double ***E1site,***E2site;
double fitcst,thrsld;
double Vspec,mspec;
char densread[50];
char iontype[50];
double part_funct,part_funct1,part_funct2;

void read();
void read_planar();
void read_porous();
void read_coarse();
void spin_echo();
void propagation();
void random_lattice();
void fourier_transform();
void autocorrel();
void equilibration();
double generlognorm(double mean, double std);
double genergauss(double mean, double std);
double energybarrier(double porsize, double temp, double lattpar, int d);
void neighbouring();
void ring_stats();
void find_rings(int length);
void find_normals();
double calc_angle(double coorda1, double coorda2, double coorda3, double coordb1, double coordb2, double coordb3);
int find_anglezone(double angle);

int **imatrix(int nl, int nc);
int ***tdimatrix (int X_SIZE, int Y_SIZE, int Z_SIZE);
double *dvector(int n);
double **dmatrix(int nl, int nc);
int *ivector(int n);
double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE);
double ****fddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE, int W_SIZE);
double *****fiveddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE, int W_SIZE, int V_SIZE);


int main(void)
{
 clock_t begin,end;
 int niter,ndone;
 double time_spent;

 printf("Hello!\n");
 
 begin=clock();

 /* Input data common to all files */
 read();
 /* Input data by simulation type */
 if(strcmp(simtype,"planar")==0) read_planar();
 if(strcmp(simtype,"porous")==0) read_porous();
 if(strcmp(simtype,"coarse")==0) read_coarse();

 tau=0;

 /* Equilibration period */
 if(strcmp(simtype,"coarse")==0)
   {
    equilibration();
   }

 /* Main function */
 propagation();

 /* Spin echo sequence */
 /* THIS SHOULD NOT BE USED AS IT HAS NOT BEEN TESTED */
 if(strcmp(pseq,"echo")==0) 
   {
    fprintf(outspinecho,"%lf	%lf\n",tau*dtime,int_tau0);
    niter=(Nvalues/2)/10;
    ndone=1;
    for(tau=0;tau<=Nvalues/2;tau+=500)
       {
        tau=Nvalues/2;
        spin_echo();
        printf("%d done on %d, tau = %lf\n",ndone,niter,tau*dtime);
        ndone++;
        fprintf(outspinecho,"%lf	%lf\n",tau*dtime,int_tau);
       }
    fclose(outspinecho);
   }

 end=clock();
 time_spent=(end-begin)/CLOCKS_PER_SEC;
 printf("Execution time: %lf.\n",time_spent);

 return 0;
}

/* Reads "non specific" input data and allocates matrices and vectors */
void read()
{
 FILE *in,*out;
 char ligne[L];
 int i,j,k,l;
 double dx,dy,dz,dist;
 double xmin,ymin,zmin;

 dim=3;							/* Lattice in three dimensions */
 dir=2*dim;
 printf("%dD lattice, %d directions\n",dim,dir); 

 scanf("%d %d %d",&Nx,&Ny,&Nz);				/* Number of sites in x and y directions */
 scanf("%lf",&a);					/* Lattice parameter */
 scanf("%lf",&dwell);					/* Dwell time */
 scanf("%d",&nsample);					/* Number of time steps between measures */
 /* The dwell time and sampling define the timestep */
 dtime=dwell/(nsample*1.0);
 /* Velocity */
 v0=a/dtime;
 /* Spectral width */
 printf("Spectral width (Hz): %lf     %lf\n",-1.0/(2.0*dwell),1.0/(2.0*dwell));
 scanf("%lf",&T);					/* Temperature of the simulation */
 scanf("%lf",&larmorfreq);				/* Larmor frequency of the studied nucleus under the magnetic field considered */
 scanf("%d",&seed);					/* Seed for the steps using random */
 /* This is for the generation of pseudo random numbers */
 srand(seed);
 /* Last parameters for the initial densities and obstacles */
 scanf("%d",&test_latt);				/* Obstacles? Lattice given? */
 if(test_latt==1) scanf("%s",namelatt);			/* Name of the file with the lattice */
 if(test_latt==2) scanf("%lf",&coverage); 		/* Fraction of the lattice covered by obstacles */
 if(test_latt==3) 
   {
    scanf("%s %lf",namepos,&Dmin); 			/* Name of the carbon positions' file and minimum distance */
   }
 scanf("%d",&Nvalues);					/* Number of values for the Fourier transform */
 scanf("%s",pseq);					/* Pulse sequence */
 scanf("%s",simtype);					/* Simulation type */

 Nmin=Nvalues*nsample;

 Lattice=tdimatrix(Nx+1,Ny+1,Nz+1); 
 Esite=tddmatrix(Nx+1,Ny+1,Nz+1); 
 alpha=fddmatrix(Nx+1,Ny+1,Nz+1,dir);
 probv0=fiveddmatrix(Nx+1,Ny+1,Nz+1,dim,2);
 /*probvo[*][*][0][*] = x and probvo[*][*][1][*] = y and probvo[*][*][2][*] = z */
 dens=fddmatrix(Nx+1,Ny+1,Nz+1,2);
 wij=tddmatrix(Nx+1,Ny+1,Nz+1);
 ImG=fddmatrix(Nx+1,Ny+1,Nz+1,2);  ReG=fddmatrix(Nx+1,Ny+1,Nz+1,2); 
 ImGtot=dvector(Nmin+1);           ReGtot=dvector(Nmin+1);
 ImGtotsmp=dvector(Nmin+1);        ReGtotsmp=dvector(Nmin+1);
 ImGcf=dvector(Nmin+1);            ReGcf=dvector(Nmin+1);
 ImGsmpcf=dvector(Nmin+1);         ReGsmpcf=dvector(Nmin+1);
 xvacf=dvector(Nz+1);  xvacfx=dvector(Nz+1);  xvacfy=dvector(Nz+1);  xvacfz=dvector(Nz+1);
 yvacf=dvector(Nz+1);  yvacfx=dvector(Nz+1);  yvacfy=dvector(Nz+1);  yvacfz=dvector(Nz+1);
 zvacf=dvector(Nz+1);  zvacfx=dvector(Nz+1);  zvacfy=dvector(Nz+1);  zvacfz=dvector(Nz+1);
 xsumvacf=dvector(Nx+1);  xsumvacfx=dvector(Nx+1);  xsumvacfy=dvector(Nx+1);  xsumvacfz=dvector(Nx+1);
 ysumvacf=dvector(Ny+1);  ysumvacfx=dvector(Ny+1);  ysumvacfy=dvector(Ny+1);  ysumvacfz=dvector(Ny+1);
 zsumvacf=dvector(Nz+1);  zsumvacfx=dvector(Nz+1);  zsumvacfy=dvector(Nz+1);  zsumvacfz=dvector(Nz+1);
 sitesize=tddmatrix(Nx+1,Ny+1,Nz+1);
 pvol=tddmatrix(Nx+1,Ny+1,Nz+1);
 vol_occ=tddmatrix(Nx+1,Ny+1,Nz+1);
 Ntot=tddmatrix(Nx+1,Ny+1,Nz+1);
 pVtot=tddmatrix(Nx+1,Ny+1,Nz+1);
 msdvacf=dvector(Nmin+1);
 msdvacfx=dvector(Nmin+1);
 msdvacfy=dvector(Nmin+1);
 msdvacfz=dvector(Nmin+1);
 dens1=fddmatrix(Nx+1,Ny+1,Nz+1,2);
 dens2=fddmatrix(Nx+1,Ny+1,Nz+1,2);
 E1site=tddmatrix(Nx+1,Ny+1,Nz+1);
 E2site=tddmatrix(Nx+1,Ny+1,Nz+1);
 inpore_c=tddmatrix(Nx+1,Ny+1,Nz+1);
 pinpore_c=tddmatrix(Nx+1,Ny+1,Nz+1);
 surface_pore=tddmatrix(Nx+1,Ny+1,Nz+1);
 

 if(strcmp(pseq,"echo")==0) 
   {
    outspinecho=fopen("Intensity_tau.dat","w");
    fprintf(outspinecho,"# Tau - Intensity\n");
   }

 /* Nl =  total number of lattice sites 
    Nf =  total number of free sites */
 Nl=Nx*Ny*Nz;
 /* Lattice generated with no obstacles */
 if(test_latt==0) {Nf=Nl;}
 /* Lattice generated with a network of obstacles given as an input */
 if(test_latt==1) 
   {
    Nf=0;
    in=fopen(namelatt,"r");
    for(i=1;i<=Nx;i++)
       {
	for(j=1;j<=Ny;j++)
	   {
	    for(k=1;k<=Nz;k++)
	       {
		fgets(ligne,L,in);
		sscanf(ligne,"%d",&Lattice[i][j][k]);
		if(Lattice[i][j][k]==0) Nf+=1;
	       }
	   }
       }
    fclose(in);
    printf("%d free sites\n",Nf);
    printf("%d obstacles\n",Nl-Nf);
   }
 /* Lattice generated with randomly distributed obstacles */
 if(test_latt==2)
   {
    Nf=0;
    out=fopen("Network_random.dat","w");
    random_lattice();
    for(k=1;k<=Nz;k++)
       {
	for(j=1;j<=Ny;j++)
	   {
	    for(i=1;i<=Nx;i++)
	       {
		fprintf(out,"%d\n",Lattice[i][j][k]);
		if(Lattice[i][j][k]==0) Nf+=1;
	       }
	   } 
       }
    fclose(out);
    printf("%d free sites\n",Nf);
    printf("%d obstacles\n",Nl-Nf);
    printf("Obtained coverage: %lf\n",(Nl-Nf)*1.0/(Nl*1.0));
   }
 /* Lattice generated with a network of obstacles given as an input from carbon structures */
 if(test_latt==3)
   {
    Nf=0;
    in=fopen(namepos,"r");
    /* Number of atoms in the xyz file */
    fgets(ligne,L,in);
    sscanf(ligne,"%d",&natoms);
    /* Skip comment line */
    fgets(ligne,L,in);
    X=dvector(natoms+1);
    Y=dvector(natoms+1);
    Z=dvector(natoms+1);
    /* Reading all the atomic positions */
    for(i=1;i<=natoms;i++)
       {
	fgets(ligne,L,in);
	sscanf(ligne,"%*s %lf %lf %lf",&X[i],&Y[i],&Z[i]);
	if(i==1) {xmin=X[i]; ymin=Y[i]; zmin=Z[i];}
	if(i>1) 
	  {
	   if(xmin>=X[i]) xmin=X[i];
	   if(ymin>=Y[i]) ymin=Y[i];
	   if(zmin>=Z[i]) zmin=Z[i];
	  }
       }
    fclose(in);
    out=fopen(namepos,"w");
    fprintf(out,"%d\n",natoms);
    fprintf(out,"Modified xyz file\n");
    /* Shift the coordinates so that all coordinates are positive */
    for(i=1;i<=natoms;i++)
       {
	X[i]-=xmin;
	Y[i]-=ymin;
	Z[i]-=zmin;
    	fprintf(out,"C	%lf	%lf	%lf\n",X[i],Y[i],Z[i]);
       }
    fclose(out);
    Lx=a*Nx;
    Ly=a*Ny;
    Lz=a*Nz;
    out=fopen("Network_acc.dat","w");
    /* Check which sites are accessible and which are not */
    for(i=1;i<=Nx;i++)
       {
	for(j=1;j<=Ny;j++)
	   {
	    for(k=1;k<=Nz;k++)
	       {
		l=1;
                while((Lattice[i][j][k]==0)&&(l<natoms))
                   {
                    dx=(i-0.5)*a-X[l];
                    dy=(j-0.5)*a-Y[l];
                    dz=(k-0.5)*a-Z[l];
                    /* Periodic boundary conditions */
                    dx/=Lx;
                    dy/=Ly;
                    dz/=Lz;
                    dx-=round(dx);
                    dy-=round(dy);
                    dz-=round(dz);
                    dx*=Lx;
                    dy*=Ly;
                    dz*=Lz;
                    dist=dx*dx+dy*dy+dz*dz;
                    dist=sqrt(dist);
                    if(dist<Dmin) Lattice[i][j][k]=1;
                    l++;
                   }
		fprintf(out,"%d\n",Lattice[i][j][k]);
		if(Lattice[i][j][k]==0) Nf+=1;
	       }
	   }
       }
    fclose(out);
    printf("%d free sites\n",Nf);
    printf("%d obstacles\n",Nl-Nf);
    printf("Obtained coverage: %lf\n",(Nl-Nf)*1.0/(Nl*1.0));
   }

}


/* Reads input data for "planar" simulation type */
void read_planar()
{
 FILE *in,*out,*outxyz;
 char ligne[L];
 int i,j,k,kbis,l,ener,idist,zcount;
 double denstot,part_func;
 double rmax,dist,distmin,meanengz,meanfreqz,meandensz;

 denstot=0.0;
 part_func=0.0;

 scanf("%s",nameeng);					/* Name of the file with free energies*/
 scanf("%d",&Neng);					/* Number of different free energies */
 in=fopen(nameeng,"r");
 free_eng=dmatrix(Neng+1,2);
 for(i=1;i<=Neng;i++) 
    {
     fgets(ligne,L,in);
     sscanf(ligne,"%lf %lf",&free_eng[i][0],&free_eng[i][1]); 
    } 
 fclose(in);
 scanf("%s",poreshift);					/* Direct shifts or pore calculation or Xing model (see below)? */
 if((strcmp(poreshift,"direct")==0)||(strcmp(poreshift,"indirect")==0))
   {
    scanf("%s",namefreq1);				/* Name of the file with frequencies */
    scanf("%d",&Nfreq);					/* Number of different frequencies */
    in=fopen(namefreq1,"r");
    freq=dmatrix(Nfreq+1,2);
    for(i=1;i<=Nfreq;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf",&freq[i][0],&freq[i][1]);  
       } 
    fclose(in);
   }

 outxyz=fopen("Lattice.xyz","w");
 fprintf(outxyz,"%d\n",Nx*Ny*Nz);
 fprintf(outxyz,"Step 1\n");
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==1) fprintf(outxyz,"O	%lf	%lf	%lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
	     if(Lattice[i][j][k]==0)
	      {
	       fprintf(outxyz,"F	%lf	%lf	%lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
	       /* energy as a function of the distance to the surface */
	       idist=1;
	       distmin=fabs(free_eng[1][0]-a*(k-1));
	       for(l=2;l<=Neng;l++)
		  {
		   if(k<=Nz/2) dist=fabs(free_eng[l][0]-a*(k-1));
		   if(k>Nz/2) dist=fabs(free_eng[l][0]-a*(Nz-k));
		   if(dist<=distmin) {distmin=dist; idist=l;}
		  }
	       Esite[i][j][k]=free_eng[idist][1];
	       part_func+=exp(-Esite[i][j][k]/(kB*T));
	       /* frequency as a function of the distance to the surface */
               if(strcmp(poreshift,"direct")==0)
	         {
	          idist=1; 
	          distmin=fabs(freq[1][0]-a*(k-1));
		  for(l=2;l<=Nfreq;l++)
		     {
		      if(k<=Nz/2) dist=fabs(freq[l][0]-a*(k-1));
		      if(k>Nz/2) dist=fabs(freq[l][0]-a*(Nz-k));
		      if(dist<=distmin) {distmin=dist; idist=l;}
		     }
                  wij[i][j][k]=freq[idist][1];
	 	 }
               if(strcmp(poreshift,"indirect")==0)
	         {
		  /* Shift due to left pore wall */
	          idist=1; 
	          distmin=fabs(freq[1][0]-a*(k-1));
		  for(l=2;l<=Nfreq;l++)
		     {
		      dist=fabs(freq[l][0]-a*(k-1));
		      if(dist<=distmin) {distmin=dist; idist=l;}
		     }
                  wij[i][j][k]=freq[idist][1];
		  /* Shift due to right pore wall */
	          idist=1; 
		  kbis=Nz-k;
	          distmin=fabs(freq[1][0]-a*kbis);
		  for(l=2;l<=Nfreq;l++)
		     {
		      dist=fabs(freq[l][0]-a*kbis);
		      if(dist<=distmin) {distmin=dist; idist=l;}
		     }
                  wij[i][j][k]+=freq[idist][1];
	         }
	       /* Fit used in Xing14.pdf: Carbon, 77, 1132 (2014). */
	       if(strcmp(poreshift,"xing")==0)
		 {
		  if(k<=Nz/2) dist=a*(k-1);
		  if(k>Nz/2) dist=a*(Nz-k);
		  wij[i][j][k]=-24.5*exp(-1.0*pow(dist/2.27,0.754));
		 }
               /* From ppm to Hz */
               wij[i][j][k]=(wij[i][j][k]*larmorfreq);
               /* From Hz to rad.s-1 */
	       wij[i][j][k]*=2.0*PI;
	     }
          }
      }
  }
 fclose(outxyz); 
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
	        /* Initialisation of the density to Boltzmann densities and 0.0 for obstacles */
	        dens[i][j][k][0]=Nf*exp(-Esite[i][j][k]/(kB*T))/(part_func);
	        denstot+=dens[i][j][k][0];
	       }
	    }
	}
   }

 out=fopen("Profiles.dat","w"); 
 fprintf(out,"# z - Energy - Frequency - Density\n");
 for(k=1;k<=Nz;k++)
    {
     meanengz=0.0;
     meanfreqz=0.0;
     meandensz=0.0;
     zcount=0.0;
     for(j=1;j<=Ny;j++)
	{
	 for(i=1;i<=Nx;i++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
		meanengz+=Esite[i][j][k];
	        meanfreqz+=((wij[i][j][k]/2.0*PI))/larmorfreq;
	  	meandensz+=dens[i][j][k][0];
		zcount++; 
	       } 
	    }
	}
     if(zcount!=0) {fprintf(out,"%e	%lf	%lf	%lf\n",(k-0.5)*a,meanengz/(zcount*1.0),meanfreqz/(zcount*1.0),meandensz/(zcount*1.0));}
   }
 
 fclose(out);
 
 printf("Reading OK!\n");

}


/* Reads input data for "porous" simulation type */
void read_porous()
{
 FILE *in,*out,*out2,*outxyz;
 char ligne[L];
 int i,j,k,n,ibis,jbis,kbis,l,ener,idist,zcount,ibin,lsite;
 double denstot,part_func,dx,dy,dz,Lx,Ly,Lz,rmax,dist,distmin,distsurf;
 double meanengz,meanfreqz,meandensz,dr,distcarbon,distsite,Dxx,Dyy,Dzz,normn,costheta;
 double *distdenscarbon,*Xeng,*Yeng,*Zeng,*readeng,*Xshift,*Yshift,*Zshift,*readshift;

 denstot=0.0;
 part_func=0.0;

 /* Preparation for the calculation of density as a function of distance to the carbon / carbon surface */
 dr=a/2.0; 
 distdenscarbon=dvector(Nz+1);
 for(k=1;k<=Nz;k++) {distdenscarbon[k]=0.0;}

 Lx=a*Nx;
 Ly=a*Ny;
 Lz=a*Nz;

 scanf("%s",test_feng);					/* Type of free energies assignement */
 scanf("%s %d",nameeng,&Neng);				/* Name of the free energies file and number of values */
 scanf("%s",test_shift);				/* Type of chemical shifts assignement */
 if((strcmp(test_shift,"1D")==0)||(strcmp(test_shift,"3D")==0))
   {
    scanf("%s %d",namefreq1,&Nfreq);			/* Name of the chemical shifts file and number of values */
   }

 /* Assignment of free energies depending on model */
 /* Free energies from a 1D eng = f(distance) model */
 if(strcmp(test_feng,"1D")==0)
   {
    /* Read free energies from file */
    in=fopen(nameeng,"r");
    free_eng=dmatrix(Neng+1,2);
    for(i=1;i<=Neng;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf",&free_eng[i][0],&free_eng[i][1]); 
       } 
    fclose(in);
    outxyz=fopen("Lattice.xyz","w");
    fprintf(outxyz,"%d\n",Nx*Ny*Nz);
    fprintf(outxyz,"Step 1\n");
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==1) fprintf(outxyz,"O  %lf     %lf     %lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
                if(Lattice[i][j][k]==0)
                  {
               	   fprintf(outxyz,"F        %lf     %lf     %lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
                   /* Loop over carbon atoms to find the closest */
                   distcarbon=Nz*dr;
                   for(l=1;l<=natoms;l++)
                      {
                       dx=(i-0.5)*a-X[l]; dy=(j-0.5)*a-Y[l]; dz=(k-0.5)*a-Z[l];
                       /* Periodic boundary conditions */
                       dx/=Lx; dy/=Ly; dz/=Lz;
                       dx-=round(dx); dy-=round(dy); dz-=round(dz);
                       dx*=Lx; dy*=Ly; dz*=Lz;
                       dist=dx*dx+dy*dy+dz*dz;
                       dist=sqrt(dist);
                       if(dist<distcarbon) distcarbon=dist;
                      }
                  /* energy as a function of the distance to the carbon */
                  idist=1;
                  distmin=fabs(free_eng[1][0]-distcarbon);
                  for(l=2;l<=Neng;l++)
                     {
                      dist=fabs(free_eng[l][0]-distcarbon);
                      if(dist<=distmin) {distmin=dist; idist=l;}
                     }
                  Esite[i][j][k]=free_eng[idist][1];
                  part_func+=exp(-Esite[i][j][k]/(kB*T));
                  /* Density as a function of distance to the closest carbon */
                  ibin=floor(distcarbon/dr);
                  if((ibin>0)&&(ibin<=Nz))
                    {
                     distdenscarbon[ibin]+=exp(-Esite[i][j][k]/(kB*T));
                    }
	         }
	      }
           }
        }
   fclose(outxyz); 
  }

 /* Free energies from a 3D matrix */
 if(strcmp(test_feng,"3D")==0)
   {
    Xeng=dvector(Neng+1);
    Yeng=dvector(Neng+1);
    Zeng=dvector(Neng+1);
    readeng=dvector(Neng+1);
    /* Read free energies from file */
    in=fopen(nameeng,"r");
    for(i=1;i<=Neng;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf %lf %lf",&Xeng[i],&Yeng[i],&Zeng[i],&readeng[i]); 
       } 
    fclose(in);
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
                   /* Loop over read site to find the closest */
                   distsite=Nz*dr;
		   lsite=0;
                   for(l=1;l<=Neng;l++)
                      {
                       dx=(i-0.5)*a-Xeng[l]; dy=(j-0.5)*a-Yeng[l]; dz=(k-0.5)*a-Zeng[l];
                       /* Periodic boundary conditions */
                       dx/=Lx; dy/=Ly; dz/=Lz;
                       dx-=round(dx); dy-=round(dy); dz-=round(dz);
                       dx*=Lx; dy*=Ly; dz*=Lz;
                       dist=dx*dx+dy*dy+dz*dz;
                       dist=sqrt(dist);
                       if(dist<distsite) {distsite=dist; lsite=l;}
		      }
                    Esite[i][j][k]=readeng[lsite];
                    part_func+=exp(-Esite[i][j][k]/(kB*T));
		  }
	       }
	   }
       }
   }

	  
 /* Assignment of chemical shifts depending on model */
 /* Shifts from a 1D eng = f(distance) model */
 if(strcmp(test_shift,"1D")==0)
   {
    /* Read chemical shifts from file */
    in=fopen(namefreq1,"r");
    freq=dmatrix(Nfreq+1,2);
    for(i=1;i<=Nfreq;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf",&freq[i][0],&freq[i][1]);  
       } 
    fclose(in);
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
                   /* Loop over carbon atoms to find the closest */
                   distcarbon=Nz*dr;
                   for(l=1;l<=natoms;l++)
                      {
                       dx=(i-0.5)*a-X[l]; dy=(j-0.5)*a-Y[l]; dz=(k-0.5)*a-Z[l];
                       /* Periodic boundary conditions */
                       dx/=Lx; dy/=Ly; dz/=Lz;
                       dx-=round(dx); dy-=round(dy); dz-=round(dz);
                       dx*=Lx; dy*=Ly; dz*=Lz;
                       dist=dx*dx+dy*dy+dz*dz;
                       dist=sqrt(dist);
                       if(dist<distcarbon) distcarbon=dist;
                      }
	           /* frequency as a function of the distance to the surface */
                   if(strcmp(test_shift,"1D")==0)
	             {
	              idist=1;
	              distmin=fabs(freq[1][0]-distcarbon);
	              for(l=2;l<=Nfreq;l++)
		         {
		          dist=fabs(freq[l][0]-distcarbon);
		          if(dist<=distmin) {distmin=dist; idist=l;}
		         }
                      wij[i][j][k]=freq[idist][1];
                      /* From ppm to Hz */
                      wij[i][j][k]=(wij[i][j][k]*larmorfreq);
                      /* From Hz to rad.s-1 */
	              wij[i][j][k]*=2.0*PI;
	 	     }
	         }
	      }
           }
        }
  }

 /* Assignment of chemical shifts depending on model */
 /* Shifts from the Xing model */
 if(strcmp(test_shift,"Xing")==0)
   {
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
                   /* Loop over carbon atoms to find the closest */
                   distcarbon=Nz*dr;
                   for(l=1;l<=natoms;l++)
                      {
                       dx=(i-0.5)*a-X[l]; dy=(j-0.5)*a-Y[l]; dz=(k-0.5)*a-Z[l];
                       /* Periodic boundary conditions */
                       dx/=Lx; dy/=Ly; dz/=Lz;
                       dx-=round(dx); dy-=round(dy); dz-=round(dz);
                       dx*=Lx; dy*=Ly; dz*=Lz;
                       dist=dx*dx+dy*dy+dz*dz;
                       dist=sqrt(dist);
                       if(dist<distcarbon) distcarbon=dist;
                      }
	           /* frequency as a function of the distance to the surface */
	           /* Fit used in Xing14.pdf: Carbon, 77, 1132 (2014). */
		   wij[i][j][k]=-24.5*exp(-1.0*pow(distsurf/2.27,0.754));
                   /* From ppm to Hz */
                   wij[i][j][k]=(wij[i][j][k]*larmorfreq);
                   /* From Hz to rad.s-1 */
	           wij[i][j][k]*=2.0*PI;
	         }
	      }
           }
        }
  }

 /* Shifts from a 3D matrix */
 if(strcmp(test_shift,"3D")==0)
   {
    Xshift=dvector(Nfreq+1);
    Yshift=dvector(Nfreq+1);
    Zshift=dvector(Nfreq+1);
    readshift=dvector(Nfreq+1);
    /* Read free energies from file */
    in=fopen(namefreq1,"r");
    for(i=1;i<=Nfreq;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf %lf %lf",&Xshift[i],&Yshift[i],&Zshift[i],&readshift[i]); 
       } 
    fclose(in);
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
                   /* Loop over read site to find the closest */
                   distsite=Nz*dr;
		   lsite=0;
                   for(l=1;l<=Nfreq;l++)
                      {
                       dx=(i-0.5)*a-Xshift[l]; dy=(j-0.5)*a-Yshift[l]; dz=(k-0.5)*a-Zshift[l];
                       /* Periodic boundary conditions */
                       dx/=Lx; dy/=Ly; dz/=Lz;
                       dx-=round(dx); dy-=round(dy); dz-=round(dz);
                       dx*=Lx; dy*=Ly; dz*=Lz;
                       dist=dx*dx+dy*dy+dz*dz;
                       dist=sqrt(dist);
                       if(dist<distsite) {distsite=dist; lsite=l;}
		      }
                   wij[i][j][k]=readshift[lsite];
                   /* From ppm to Hz */
                   wij[i][j][k]=(wij[i][j][k]*larmorfreq);
                   /* From Hz to rad.s-1 */
	           wij[i][j][k]*=2.0*PI;
		  }
	       }
	   }
       }
   }

 /* Shifts from the dipolar model */
 if(strcmp(test_shift,"dipolar")==0)
   {
    /* maximum length of the rings in number of atoms */
    Neigh=imatrix(natoms+1,natoms+1);
    Neigh_modif=imatrix(natoms+1,natoms+1);
    current_path=ivector(Lmax+1);
    Histo_rings=ivector(Lmax+1);
    PATHS=imatrix(3*natoms+1,Lmax+1);
    Histo_dist=ivector(nbins+1);
    /* For the distances histogram */
    dmax=Lx;
    if(Ly>=dmax) dmax=Ly;
    if(Lz>=dmax) dmax=Lz;
    dmax=dmax*sqrt(3.0)+1.0;                       /* Longest distance = cube diagonal */
    neighbouring();
    ring_stats();
    find_normals();
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
		   Dxx=0.0; Dyy=0.0; Dzz=0.0;
		   for(n=1;n<=Nrings;n++)
		      {
		       dx=(i-0.5)*a-xcom[n];
		       dy=(j-0.5)*a-ycom[n];
		       dz=(k-0.5)*a-zcom[n];
	 	       /* Periodic boundary conditions */
	    	       if(dx>(Lx/2.0)) dx-=Lx;
            	       if(dx<(-Lx/2.0)) dx+=Lx;
	    	       if(dy>(Ly/2.0)) dy-=Ly;
            	       if(dy<(-Ly/2.0)) dy+=Ly;
	    	       if(dz>(Lz/2.0)) dz-=Lz;
            	       if(dz<(-Lz/2.0)) dz+=Lz;
		       dist=dx*dx+dy*dy+dz*dz;
		       dist=sqrt(dist);
		       if(dist<=cutoffshift)
			 {
			  normn=xnorm[n]*xnorm[n]+ynorm[n]*ynorm[n]+znorm[n]*znorm[n];
			  normn=sqrt(normn);
			  costheta=dx*xnorm[n]+dy*ynorm[n]+dz*znorm[n];
			  costheta/=(dist*normn);
			  Dxx+=Chiperp*(1.0/4.0*PI)*((3.0*sin(acos(costheta))*sin(acos(costheta))-1.0)/(dist*dist*dist));
			  Dyy-=Chiperp*(1.0/4.0*PI)*(1.0/(dist*dist*dist));
			  Dzz+=Chipara*(1.0/4.0*PI)*((3*costheta*costheta-1.0)/(dist*dist*dist));
			 }
		      }
                   wij[i][j][k]=(1.0/3.0)*(Dzz+Dyy+Dxx);
                   /* From ppm to Hz */
                   wij[i][j][k]=(wij[i][j][k]*larmorfreq);
                   /* From Hz to rad.s-1 */
	           wij[i][j][k]*=2.0*PI;
		  }
	       }
	   }
       }
   }

 /* Initialisation of the density to Boltzmann densities and 0.0 for obstacles */
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
	        dens[i][j][k][0]=Nf*exp(-Esite[i][j][k]/(kB*T))/(part_func);
	        denstot+=dens[i][j][k][0];
	       }
	    }
	}
   }

 out=fopen("Profiles.dat","w"); 
 fprintf(out,"# z - Energy - Frequency\n");
 for(k=1;k<=Nz;k++)
    {
     meanengz=0.0; meanfreqz=0.0; meandensz=0.0; zcount=0.0;
     for(j=1;j<=Ny;j++)
	{
	 for(i=1;i<=Nx;i++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
		meanengz+=Esite[i][j][k]; meanfreqz+=((wij[i][j][k]/2.0*PI))/larmorfreq; meandensz+=dens[i][j][k][0]; zcount++; 
	       } 
	    }
	}
     if(zcount!=0) {fprintf(out,"%e	%lf	%lf	%lf\n",k*a,meanengz/(zcount*1.0),meanfreqz/(zcount*1.0),meandensz/(zcount*1.0));}
   }
 fclose(out);
 
 out=fopen("Density_wrt_surface.dat","w"); 
 out2=fopen("Density_wrt_carbonpos.dat","w"); 
 fprintf(out,"# dist - Density\n");
 fprintf(out2,"# dist - Density\n");
 for(i=1;i<=Nz;i++)
    {
     fprintf(out2,"%e	%lf\n",i*a,distdenscarbon[i]*Nf/(part_func)); 
    }

 fclose(out);
 fclose(out2);

 printf("Reading OK!\n");

}


/* Reads input data for "coarse" simulation type */
void read_coarse()
{
 FILE *in,*inm,*out,*outxyz,*outpores,*outbarr,*outshifts,*outpsd, *outpsdpart;
 FILE *outgradpore,*outenergy,*outvolume,*outmass;
 char ligne[L];
 char cover_pref[10],namedens_an[50],namedens_cat[50],namedens_solv[50];
 char iontype[50],solvent_activ[50],psdread[50],namepsdm[50];
 int i,j,k,l,zcount,count,countb,i1,i2,j1,j2,k1,k2,nlinespsd,nlinessurf,nlinespsdm;
 int ipore,jbefore,jafter,gradient,igrad,jgrad,kgrad,rinumber,ifind,countt,idist,jdist,kdist;
 double denstot,part_func,poresurf,intrap,slope,gradmin,gradcomp,rdnumber;
 double poresize,meanengz,meanfreqz,meandensz,distmin,sumpb,sumpbm,distcentpart;
 double *gradw,*gradD,*pospsd,*pbapsd,*possurf,*pbasurf;
 double *pospsdm,*pbapsdm,poremax,porecomp,poremin,rdmix;
 double Nanion,Ncation,Nsolvent,Ntotanion,Ntotcation,Ntotsolvent;
 double sitei,sitej,Ea_i,Ea_j,min_dens,avg_poresize,sumporesize,porevol,totalvol;
 double aniondens,cationdens,solventdens,monodens;
 double anionmass,cationmass,solventmass;
 double totanmass,totcatmass,totsolvmass;
 double poremass_anion,poremass_cation,poremass_solvent,poremass;
 double avg_intdens,elecsurf,volt,cpt,charge,rhotrshld,Nsolvat,mixratio;
 double *gradp,*gradN,*gradD1,*gradD2,*gradVpore,*gradVocc,*gradcc,*gradpcc,*fillset;

 part_func=0.0;
 part_funct1=0.0;
 part_funct2=0.0;

 scanf("%d",&Npores);                                      /* Number of different pore sizes in the input files */
 scanf("%s",densread);                                     /* can take distinguished or undistinguished values */
 if(strcmp(densread,"distinguished")==0)
   {
     scanf("%s",namedens_an);                               /* Name of the file with anion integrated densities */
     scanf("%s",namedens_cat);                              /* Name of the file with cation integrated densities */
     scanf("%s",solvent_activ);                             /* Presence (solvent_ON) or abscence (solvent_OFF) of solvent */
     if(strcmp(solvent_activ,"solvent_ON")==0)
       {
        scanf("%s",namedens_solv);                             /* Name of the file with solvent integrated densities */
       }
     intdens=dmatrix(Npores+1,4);
     in=fopen(namedens_an,"r");
     for(i=1;i<=Npores;i++)
        {
         fgets(ligne,L,in);
         sscanf(ligne,"%lf %lf",&intdens[i][0],&intdens[i][1]);
        }
     fclose(in);
     in=fopen(namedens_cat,"r");
     for(i=1;i<=Npores;i++)
        {
         fgets(ligne,L,in);
         sscanf(ligne,"%*lf %lf",&intdens[i][2]);
        }
     fclose(in);
     if(strcmp(solvent_activ,"solvent_ON")==0)
       {
        in=fopen(namedens_solv,"r");
        for(i=1;i<=Npores;i++)
           {
            fgets(ligne,L,in);
            sscanf(ligne,"%*lf %lf",&intdens[i][3]);
           }
        fclose(in);
       }
     if(strcmp(solvent_activ,"solvent_OFF")==0)
       {
        for(i=1;i<=Npores;i++)
           {
            intdens[i][3]=0.0;
           }
       }
     scanf("%s",iontype);                                    /* Type of mol (= cation, anion or solvent) for which the calculation will be performed */ 
     if(strcmp(solvent_activ,"solvent_OFF")==0&&strcmp(iontype,"solvent")==0)
       {
        printf("ERROR: Cannot choose solvent - solvent is not considered\n");
        exit(0);
       }
   }
 if(strcmp(densread,"undistinguished")==0)
   {
    scanf("%s",namedens);					/* Name of the file with integrated densities */
    intdens=dmatrix(Npores+1,2);
    in=fopen(namedens,"r");
    for(i=1;i<=Npores;i++)
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf",&intdens[i][0],&intdens[i][1]);
       }
    fclose(in);
   }
 area=dvector(4);
 scanf("%s %lf",namefreq1,&area[1]);			/* Name of the file with frequencies and area */
 scanf("%s %lf",namefreq2,&area[2]);			/* Name of the file with frequencies and area */
 scanf("%s %lf",namefreq3,&area[3]);			/* Name of the file with frequencies and area */
 
 /* Reading of the chemical shifts (in ppm) as a function of pore size */
 freq=dmatrix(Npores+1,4);
 in=fopen(namefreq1,"r");
 for(i=1;i<=Npores;i++) 
    {
     fgets(ligne,L,in);
     sscanf(ligne,"%lf %lf",&freq[i][0],&freq[i][1]);  
    } 
 fclose(in);
 in=fopen(namefreq2,"r");
 for(i=1;i<=Npores;i++) 
    {
     fgets(ligne,L,in);
     sscanf(ligne,"%*lf %lf",&freq[i][2]);  
    } 
 fclose(in);
 in=fopen(namefreq3,"r");
 for(i=1;i<=Npores;i++) 
    {
     fgets(ligne,L,in);
     sscanf(ligne,"%*lf %lf",&freq[i][3]);  
    } 
 fclose(in);
 /* Parameters for the log normal distributions of pore sizes and pore areas */
 scanf("%s",psdtype);
 if(strcmp(psdtype,"lognormal")==0) scanf("%lf %lf",&meanpsd,&stdpsd);
 if(strcmp(psdtype,"discrete")==0) 
   {
    scanf("%s",psdread);
    if(strcmp(psdread,"mono")==0)
      {
       scanf("%s %d",namepsd,&nlinespsd);
       sumpb=0;
       pospsd=dvector(nlinespsd+1);
       pbapsd=dvector(nlinespsd+1);
       in=fopen(namepsd,"r");
       for(i=1;i<=nlinespsd;i++) 
          {
           fgets(ligne,L,in);
           sscanf(ligne,"%lf %lf",&pospsd[i],&pbapsd[i]);  
           sumpb+=pbapsd[i];
          }	
       fclose(in);
       for(i=1;i<=nlinespsd;i++) pbapsd[i]/=sumpb;
      }
    if(strcmp(psdread,"mixed")==0)
      {
       scanf("%s %d",namepsd,&nlinespsd);
       scanf("%s %d",namepsdm,&nlinespsdm);
       scanf("%lf",&mixratio); 
       sumpb=0;
       sumpbm=0;
       pospsd=dvector(nlinespsd+1);
       pbapsd=dvector(nlinespsd+1);
       pospsdm=dvector(nlinespsdm+1);
       pbapsdm=dvector(nlinespsdm+1);
       in=fopen(namepsd,"r");
       inm=fopen(namepsdm,"r");
       for(i=1;i<=nlinespsd;i++)
          {
           fgets(ligne,L,in);
           sscanf(ligne,"%lf %lf",&pospsd[i],&pbapsd[i]);
           sumpb+=pbapsd[i];
          }
       for(i=1;i<=nlinespsdm;i++)
          {
           fgets(ligne,L,inm);
           sscanf(ligne,"%lf %lf",&pospsdm[i],&pbapsdm[i]);
           sumpbm+=pbapsdm[i];
          }
       fclose(in);
       fclose(inm);
       for(i=1;i<=nlinespsd;i++) pbapsd[i]/=sumpb;
       for(i=1;i<=nlinespsdm;i++) pbapsdm[i]/=sumpbm;
      }
   }
 scanf("%s",surftype);
 if(strcmp(surftype,"lognormal")==0) scanf("%lf %lf",&meansurf,&stdsurf);
 if(strcmp(surftype,"discrete")==0) 
   {
    scanf("%s %d",namesurf,&nlinessurf);
    sumpb=0;
    possurf=dvector(nlinessurf+1);
    pbasurf=dvector(nlinessurf+1);
    in=fopen(namesurf,"r");
    for(i=1;i<=nlinessurf;i++) 
       {
        fgets(ligne,L,in);
        sscanf(ligne,"%lf %lf",&possurf[i],&pbasurf[i]);  
        sumpb+=pbasurf[i];
       }	
    fclose(in);
    for(i=1;i<=nlinessurf;i++) pbasurf[i]/=sumpb;
   }


 /* Reading the mode of calculation of the activation energies (0 for random, 1 for energies from fitting experimental data) */
 scanf("%d",&Enbar_mode); 

 /* Parameters for the Gaussian distribution of barrier heights */
 if(Enbar_mode==0) scanf("%lf %lf",&meanEa,&stdEa);

 if(Enbar_mode==1) scanf("%lf",&T_ref); /* Reference temperature at which the energy barriers are calculated */

 scanf("%d",&gradient);			/* Calculation with a gradient of pores? (0 for no, 1 for yes) */

 if(strcmp(densread,"distinguished")==0)
   {
    if(strcmp(solvent_activ,"solvent_ON")==0)
      {
       scanf("%lf %lf %lf %lf",&Vanion,&Vcation,&Vsolvent,&volt);    /* Volume of the considered anion and cation and solvent in ang3 and the cell voltage in V */
       scanf("%lf %lf %lf",&anionmass,&cationmass,&solventmass);     /* mass of the anion, cation and solvent molecules in u unit */
      }
    if(strcmp(solvent_activ,"solvent_OFF")==0)
      {
       scanf("%lf %lf %lf",&Vanion,&Vcation,&volt);
       Vsolvent=0.0;
       scanf("%lf %lf",&anionmass,&cationmass);
       solventmass=0.0;
      }
   }
 
 if(strcmp(densread,"undistinguished")==0)
   {
    scanf("%lf %lf",&Vspec,&volt);   /* Volume of the considered particule in ang3 and cell voltage in V */
    scanf("%lf",&mspec);             /* Mass of the considered particule */
   }

 scanf("%d %lf",&step,&thrsld);         /* Number of steps and threshold for the equilibration process */
 scanf("%s %d",wrtdns,&nsmple);         /* Writing of density files and number of steps between checks*/
 scanf("%lf",&rhotrshld);               /* Minimum in-pore density to be considered for calculations */  
 scanf("%s",cover_pref);                /* Putting obstacles randomly (=none) or on small pores (=small) or on large pores (=large) */ 
 scanf("%lf",&particlesize);		/* Particle size in lattice units */
 if(particlesize==0.0) printf("Warning: The lattice contains only bulk!\n");
 if((0.0<particlesize)&&(particlesize<Nx)) printf("Warning: The lattice contains both particle and bulk!\n");
 if((particlesize>=Nx)||(particlesize==-1)) printf("Warning: The lattice contains only particle!\n");

 outxyz=fopen("Lattice.xyz","w");
 fprintf(outxyz,"%d\n",Nx*Ny*Nz);
 fprintf(outxyz,"Step 1\n");
 outpores=fopen("Poresize_properties.dat","w");
 outbarr=fopen("Barrier_heights.dat","w");
 outshifts=fopen("wij.dat","w");
 outpsd=fopen("psd.dat","w");
 outpsdpart=fopen("psd_part.dat","w");

 if(strcmp(densread,"distinguished")==0)
   {
    outvolume=fopen("Volumes.dat","w");
    outmass=fopen("Masses.dat","w");
    fprintf(outvolume,"Poresize		Nanion		Ncation		Nsolvent	Volume_max	Volume_occ	max-occ\n");
    fprintf(outvolume,"\n");
    fprintf(outmass,"m_anion		m_cation	m_solvent	m_tot\n");
    fprintf(outmass,"\n");
    fillset=dvector(4);
    fillset[1]=Vanion;
    fillset[2]=Vcation;
    fillset[3]=Vsolvent;
    Ntotanion=0.0;
    Ntotcation=0.0;
    Ntotsolvent=0.0;
    totanmass=0.0;
    totcatmass=0.0;
    totsolvmass=0.0;
    elecsurf=0.0;
   }
 
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
        {
         for(k=1;k<=Nz;k++)
            {
             Ntot[i][j][k]=0.0;
            }
        }
    }

 sumporesize=0.0;
  for(i=1;i<=nlinespsd;i++)
    {
     sumporesize+=pospsd[i];
    }
 avg_poresize=sumporesize/nlinespsd;

 count=0;
 countb=0;
 poremax=0.0;
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==1) fprintf(outxyz,"O	%lf	%lf	%lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
	     if(Lattice[i][j][k]==0)
	       {
	        fprintf(outxyz,"F	%lf	%lf	%lf\n",(i-0.5)*a,(j-0.5)*a,(k-0.5)*a);
                retake:
                if(strcmp(psdtype,"lognormal")==0) 
	          {
		   /* Choose randomly a pore size on a log normal distribution */
	           /* ATTENTION, for usual input files the pore size falls approximately in the range from 4 to 15 nms */
	           /* MULTIPLICATION by 10 to go to to anstroms */
	           poresize=generlognorm(meanpsd,stdpsd)*10.0; 
                   sitesize[i][j][k]=poresize;

	           /* If out of the range given in the input files, get another number */
	           /* The mean and std should be chosen so that it does not happen often */
	           while((poresize>intdens[Npores][0])||(poresize<intdens[1][0]))
		        {
		         poresize=generlognorm(meanpsd,stdpsd)*10.0;
                         sitesize[i][j][k]=poresize;
		        }
	          }
                if(strcmp(psdtype,"discrete")==0) 
	          {
                   if(strcmp(psdread,"mono")==0)
                     {
                      rechoose: 
		      sumpb=0.0;
		      /* Choose randomly a pore size in the discrete distribution */
                      rdnumber=rand()%10000+1;
		      rdnumber=rdnumber/(10000*1.0);
		      /* Find the corresponding interval */
		      rinumber=1;
		      sumpb=pbapsd[1];
		      while(rdnumber>sumpb) {rinumber++; sumpb+=pbapsd[rinumber];} 
		      /* Choose the pore size */
                      poresize=pospsd[rinumber];
                      min_dens=intdens[1][0];
                      for(ifind=2;ifind<=Npores;ifind++)
                         {
                          if(intdens[ifind][0]<=min_dens) min_dens=intdens[ifind][0];
                         }
                      if(poresize<0.000001/*min_dens*/) goto rechoose;
                      if(test_latt==2&&strcmp(cover_pref,"small")==0) {if(poresize<avg_poresize*coverage) goto rechoose;}
                      if(test_latt==2&&strcmp(cover_pref,"large")==0) {if(poresize>=avg_poresize/coverage) goto rechoose;}
                      sitesize[i][j][k]=poresize;
                     }
                   if(strcmp(psdread,"mixed")==0)
                     {
                      rdmix=rand()%10000;
                      rdmix=rdmix/(10000*1.0);
                      if(rdmix<mixratio)
                        {
                         rechoose1:
                         sumpb=0.0;
                         rdnumber=rand()%10000+1;
                         rdnumber=rdnumber/(10000*1.0);
                         rinumber=1;
                         sumpb=pbapsd[1];
                         while(rdnumber>sumpb) {rinumber++; sumpb+=pbapsd[rinumber];}
                         poresize=pospsd[rinumber];
                         min_dens=intdens[1][0];
                         for(ifind=2;ifind<=Npores;ifind++)
                            {
                             if(intdens[ifind][0]<=min_dens) min_dens=intdens[ifind][0];
                            }
                         if(poresize<0.000001/*min_dens*/) goto rechoose1;
                         if(test_latt==2&&strcmp(cover_pref,"small")==0) {if(poresize<avg_poresize*coverage) goto rechoose1;}
                         if(test_latt==2&&strcmp(cover_pref,"large")==0) {if(poresize>=avg_poresize/coverage) goto rechoose1;}
                         sitesize[i][j][k]=poresize;
                        }
                      else
                        {
                         rechoose2:
                         sumpb=0.0;
                         rdnumber=rand()%10000+1;
                         rdnumber=rdnumber/(10000*1.0);
                         rinumber=1;
                         sumpb=pbapsdm[1];
                         while(rdnumber>sumpb) {rinumber++; sumpb+=pbapsdm[rinumber];}
                         poresize=pospsdm[rinumber];
                         min_dens=intdens[1][0];
                         for(ifind=2;ifind<=Npores;ifind++)
                            {
                             if(intdens[ifind][0]<=min_dens) min_dens=intdens[ifind][0];
                            }
                         if(poresize<0.000001/*min_dens*/) goto rechoose2;
                         if(test_latt==2&&strcmp(cover_pref,"small")==0) {if(poresize<avg_poresize*coverage) goto rechoose2;}
                         if(test_latt==2&&strcmp(cover_pref,"large")==0) {if(poresize>=avg_poresize/coverage) goto rechoose2;}
                         sitesize[i][j][k]=poresize;                         
                        }
                     }
	          }
	        /* Find the closest value of poresize in the input file */
	        distmin=fabs(poresize-intdens[1][0]);
		ipore=1;
		for(l=2;l<=Npores;l++)
		   {
		    if(fabs(poresize-intdens[l][0])<=distmin) {distmin=fabs(poresize-intdens[l][0]); ipore=l;}
		   }
                if(strcmp(surftype,"lognormal")==0) 
	          {
	           /* Choose a pore surface on a log normal distribution */
	           /* ATTENTION, for usual input files the pore surface is in nm2 */
	           /* MULTIPLICATION by 100 to go to to anstroms^2 */
	           poresurf=generlognorm(meansurf,stdsurf)*100.0; 
	           /* If out of the range given in the input files, get another number */
	           /* The mean and std should be chosen so that it does not happen often */
		   while((poresurf>area[3])||(poresurf<area[1]))
		        {
		         poresurf=generlognorm(meansurf,stdsurf)*100.0;
		        }
	          }
                if(strcmp(surftype,"discrete")==0) 
	          {
		   sumpb=0.0;
		   /* Choose randomly a pore size in the discrete distribution */
                   rinumber=rand()%10000+1;
		   rdnumber=rinumber/10000;
		   /* Find the corresponding interval */
		   rinumber=1;
		   sumpb=pbasurf[1];
		   while(rdnumber>sumpb) {rinumber++; sumpb+=pbasurf[rinumber];} 
		   /* Choose the pore size */
		   poresurf=possurf[rinumber];  
	          }
                surface_pore[i][j][k]=poresurf;
		/* DISTINGUISHED CASE */
                if(strcmp(densread,"distinguished")==0)
                  {
                    porevol=sitesize[i][j][k]*poresurf;
                    if(volt==0.0)
                      {
                       avg_intdens=(intdens[ipore][1]+intdens[ipore][2])/2.0;
                       intdens[ipore][1]=avg_intdens;
                       intdens[ipore][2]=avg_intdens;
                      }
                    aniondens=intdens[ipore][1]*poresurf;
                    cationdens=intdens[ipore][2]*poresurf;
                    solventdens=intdens[ipore][3]*poresurf;
                    if(strcmp(solvent_activ,"solvent_ON")==0)
                      {
                       if(aniondens<=rhotrshld||cationdens<=rhotrshld||solventdens<=rhotrshld) goto retake;
                      }
                    if(strcmp(solvent_activ,"solvent_OFF")==0)
                      {
                       if(aniondens<=rhotrshld||cationdens<=rhotrshld) goto retake;
                      }
                    Nanion=aniondens;
                    Ncation=cationdens;
                    Nsolvent=solventdens;
                    Nsolvat=Nsolvent/(Nanion+Ncation);
                    totalvol=Nanion*Vanion+Ncation*Vcation+Nsolvent*Vsolvent;
                    while(totalvol>porevol)
                         {
                          printf("Pore %d%d%d capacity exceeded: density regulation is in process...\n",i,j,k);
                          if(volt==0.0)
                            {
                             totalvol=totalvol-(fillset[1]+fillset[2]);
                             Nanion-=1.0;
                             aniondens-=1.0;
                             Nsolvent-=Nsolvat;
                             solventdens-=Nsolvat;
                             Ncation-=1.0;
                             cationdens-=1.0;
                             Nsolvent-=Nsolvat;
                             solventdens-=Nsolvat;
                             if(Nanion<=0)
                               {
                                totalvol+=fillset[1];
                                Nanion+=1.0;
                                aniondens+=1.0;
                                Nsolvent+=Nsolvat;
                                solventdens+=Nsolvat;
                               }
                             if(Ncation<=0)
                               {
                                totalvol+=fillset[2];
                                Ncation+=1.0;
                                cationdens+=1.0;
                                Nsolvent+=Nsolvat;
                                solventdens+=Nsolvat;
                               } 
                            }
                          if(volt<0.0)
                            {
                             totalvol-=fillset[1];
                             Nanion-=1.0;
                             aniondens-=1.0;
                             Nsolvent-=Nsolvat;
                             solventdens-=Nsolvat;
                             if(Nanion<=0.0)
                               {
                                totalvol+=fillset[1];
                                totalvol-=fillset[2];
                                Nanion+=1.0;
                                aniondens+=1.0;
                                Ncation-=1.0;
                                cationdens-=1.0;
                               }
                            }
                          if(volt>0.0)
                            {
                             totalvol-=fillset[2];
                             Ncation-=1.0;
                             cationdens-=1.0;
                             Nsolvent-=Nsolvat;
                             solventdens-=Nsolvat;
                             if(Ncation<=0.0)
                               {
                                totalvol+=fillset[2];
                                totalvol-=fillset[1];
                                Ncation+=1.0;
                                cationdens+=1.0;
                                Nanion-=1.0;
                                aniondens-=1.0;
                               }
                            }
                           totalvol=Nanion*Vanion+Ncation*Vcation+Nsolvent*Vsolvent;  
                          }
                    pvol[i][j][k]=porevol;
                    vol_occ[i][j][k]=totalvol;
                    inpore_c[i][j][k]=(((Nanion+Ncation)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                    fprintf(outvolume,"%lf		%lf	%lf	%lf	%lf	%lf	%lf\n",sitesize[i][j][k],Nanion,Ncation,Nsolvent,porevol,totalvol,porevol-totalvol);
                    Ntot[i][j][k]=Nanion+Ncation+Nsolvent;
                    Ntotanion+=Nanion;
                    Ntotcation+=Ncation;
                    Ntotsolvent+=Nsolvent;
                    elecsurf+=poresurf;
                    poremass_anion=Nanion*anionmass; poremass_cation=Ncation*cationmass; poremass_solvent=Nsolvent*solventmass;
                    poremass=poremass_anion+poremass_cation+poremass_solvent;
                    totanmass+=poremass_anion; totcatmass+=poremass_cation; totsolvmass+=poremass_solvent;
                    fprintf(outmass,"%lf	%lf	%lf	%lf\n",poremass_anion,poremass_cation,poremass_solvent,poremass);
                  }
		/* UNDISTINGUISHED CASE */
		if(strcmp(densread,"undistinguished")==0) 
	          {
		   monodens=intdens[ipore][1]*poresurf;
 	 	   if(monodens<=rhotrshld) goto retake;		   
		  }  
		/* Find the two values of areas that surround this pore surface */
		if(poresurf<=area[2]) {jbefore=1; jafter=2;}
		if(poresurf>area[2]) {jbefore=2; jafter=3;}
	        /* Boltzmann weight on the site is proportional to the integrated density
		  corresponding to this pore size and the surface of the pore */
		
		/* If there is a bulk region, we change the densities and frequencies */
		/***********************************/
		/* Note that the pore surface is not necessarily the same for all bulk volumes */
		/***********************************/
		if(particlesize>=0.0)
		  {
		   idist=i; jdist=j; kdist=k;
		   if(idist>(Nx/2)) idist-=Nx;
		   if(jdist>(Ny/2)) jdist-=Ny;
		   if(kdist>(Nz/2)) kdist-=Nz;
		   distcentpart=idist*idist+jdist*jdist+kdist*kdist;
		   distcentpart=sqrt(distcentpart);
		   if((strcmp(densread,"undistinguished")==0)&&(distcentpart>particlesize))
		     {
		      monodens=intdens[Npores][1]*poresurf; 
		     }
		   if((strcmp(densread,"distinguished")==0)&&(distcentpart>particlesize))
		     {
                      aniondens=intdens[Npores][1]*poresurf;
                      cationdens=intdens[Npores][2]*poresurf;
                      solventdens=intdens[Npores][3]*poresurf;
		     }
		  }
		/***********************************/
                if((strcmp(densread,"distinguished")==0)&&(strcmp(iontype,"anion")==0))
                  {
                   dens[i][j][k][0]=aniondens;
                   dens1[i][j][k][0]=cationdens;
                   dens2[i][j][k][0]=solventdens;
                   part_func+=dens[i][j][k][0];
                   part_funct1+=dens1[i][j][k][0];
                   part_funct2+=dens2[i][j][k][0];
                   part_funct=part_func;
                  }
                if((strcmp(densread,"distinguished")==0)&&(strcmp(iontype,"cation")==0)) 
                  {
                   dens[i][j][k][0]=cationdens;
                   dens1[i][j][k][0]=aniondens;
                   dens2[i][j][k][0]=solventdens;
                   part_func+=dens[i][j][k][0];
                   part_funct1+=dens1[i][j][k][0];
                   part_funct2+=dens2[i][j][k][0];
                   part_funct=part_func;
                  }
                if((strcmp(densread,"distinguished")==0)&&(strcmp(iontype,"solvent")==0)) 
                  {
                   dens[i][j][k][0]=solventdens;
                   dens1[i][j][k][0]=aniondens;
                   dens2[i][j][k][0]=cationdens;
                   part_func+=dens[i][j][k][0];
                   part_funct1+=dens1[i][j][k][0];
                   part_funct2+=dens2[i][j][k][0];
                   part_funct=part_func;
                  }
                if(strcmp(densread,"undistinguished")==0)
                  {
                   pvol[i][j][k]=sitesize[i][j][k]*poresurf; 
		   dens[i][j][k][0]=monodens;
		   Ntot[i][j][k]=dens[i][j][k][0];
                   inpore_c[i][j][k]=((Ntot[i][j][k]*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                   vol_occ[i][j][k]=Ntot[i][j][k]*Vspec;
                   part_func+=dens[i][j][k][0];
                   part_funct=part_func;
                  }
	        /* Chemical shift of this pore by intrapolation as well */
	        slope=(freq[ipore][jafter]-freq[ipore][jbefore])/(area[jafter]-area[jbefore]);
	        intrap=freq[ipore][jbefore]+slope*(poresurf-area[jbefore]);
                wij[i][j][k]=intrap;
		/* All bulk sites have a 0 frequency */
	   	if((particlesize>=0.0)&&(distcentpart>particlesize)) wij[i][j][k]=0.0;
		/* Alternate model - THIS IS A VERY SPECIFIC CASE - UNCOMMENT ONLY IF YOU KNOW WHAT THIS IS FOR */
		 /*if(((i+j+k)%2)==0) wij[i][j][k]=0.0;
		 if(((i+j+k)%2)!=0) wij[i][j][k]=-7.19;*/
		/* Write down pore sizes, pore areas, shifts and densities in case we want to plot histograms */
		count++;
		fprintf(outpores,"%d	%lf	%lf	%lf	%lf\n",count,poresize,poresurf,wij[i][j][k],dens[i][j][k][0]);
		fprintf(outshifts,"%lf\n",wij[i][j][k]);
		fprintf(outpsd,"%lf\n",poresize);
	   	if((particlesize>=0.0)&&(distcentpart<particlesize)) fprintf(outpsdpart,"%lf\n",poresize);
                /* From ppm to Hz */
                wij[i][j][k]=(wij[i][j][k]*larmorfreq);
                /* From Hz to rad.s-1 */
                wij[i][j][k]*=2.0*PI;
		if(poresize>poremax) poremax=poresize;
	       }
	    }
	}
    }

 if(strcmp(densread,"distinguished")==0) 
   {
    printf("Surface of the electrode = %e\n",elecsurf*1e-20);    /* surface in m2 */
    charge=(fabs(Ntotanion-Ntotcation)*1.60217646e-19)/(elecsurf*1e-16);     /* charge in C/cm2 */
    fprintf(outvolume,"Charge on the electrode = %e	Cell voltage = %lf\n",charge,volt);
    if(volt!=0.0)
      {
       cpt=(Ntotanion-Ntotcation)*1.60217646e-13/volt;
       cpt/=elecsurf*1e-16;
       fprintf(outvolume,"Capacitance = %lf\n",cpt);           /* Capacitance in uF/cm2 */
      }
    fprintf(outvolume,"Total number of: anions = %lf  	    cations = %lf		solvent = %lf\n",Ntotanion,Ntotcation,Ntotsolvent);
    fprintf(outvolume,"Total adsorbed population = %lf\n",Ntotanion+Ntotcation+Ntotsolvent);
    fprintf(outmass,"Total mass of: anions = %lf	cations = %lf		solvent = %lf\n",totanmass,totcatmass,totsolvmass);  /* mass in u unit */
    fprintf(outmass,"Total adsorbed mass = %lf\n",((totanmass+totcatmass+totsolvmass)*1.660540199e-18)/(elecsurf*1e-16));   /* total mass in ug/cm2 */
    fclose(outmass);
    fclose(outvolume); 
   }

 fclose(outshifts);
 if(gradient==1)
   {
    gradw=dvector(Nf+1);
    gradD=dvector(Nf+1);
    gradD1=dvector(Nf+1);
    gradD2=dvector(Nf+1);
    gradp=dvector(Nf+1);
    gradN=dvector(Nf+1);
    gradVpore=dvector(Nf+1);
    gradVocc=dvector(Nf+1);
    gradcc=dvector(Nf+1);
    gradpcc=dvector(Nf+1);
    for(i=1;i<=Nf;i++) {gradw[i]=0.0; gradD[i]=0.0; gradD1[i]=0.0; gradD2[i]=0.0; gradp[i]=0.0; gradN[i]=0.0;}
    for(i=1;i<=Nf;i++) {gradcc[i]=0.0;gradpcc[i]=0.0; gradVpore[i]=0.0; gradVocc[i]=0.0;}
    outshifts=fopen("wij.dat","w");
    outgradpore=fopen("poresize_gradient_distribution.dat","w");

    /* Sorting out the values of site properties */
    count=1;
    for(l=1;l<=Nf;l++)
       {
        poremin=poremax;
        for(i=1;i<=Nx;i++)
           {
            for(j=1;j<=Ny;j++)
               {
                for(k=1;k<=Nz;k++)
                   {
		    if(Lattice[i][j][k]==0)
		      {
		       porecomp=sitesize[i][j][k];
                       if(porecomp<=poremin)
                         {
                          igrad=i;
                          jgrad=j;
                          kgrad=k;
                          poremin=porecomp;
                         }
		      }
                   }
               }
           }
        gradw[count]=wij[igrad][jgrad][kgrad];
        gradD[count]=dens[igrad][jgrad][kgrad][0];
        if(strcmp(densread,"distinguished")==0)
          {
           gradD1[count]=dens1[igrad][jgrad][kgrad][0];
           gradD2[count]=dens2[igrad][jgrad][kgrad][0];
          }
        gradp[count]=sitesize[igrad][jgrad][kgrad];
        gradN[count]=Ntot[igrad][jgrad][kgrad];
        gradVpore[count]=pvol[igrad][jgrad][kgrad];
        gradVocc[count]=vol_occ[igrad][jgrad][kgrad];
        gradcc[count]=inpore_c[igrad][jgrad][kgrad];
        gradpcc[count]=surface_pore[igrad][jgrad][kgrad];
        count++;
        sitesize[igrad][jgrad][kgrad]=poremax*100.0;
       }
     
    /* Gradiently organising property values in case of even lattice dimensions */
    if(Nx%2==0&&Ny%2==0&&Nz%2==0)
      {
       count=1;
       countt=Nf;
       for(i=1;i<=Nx;i++)
          {
          for(j=1;j<=Ny;j++)
             {
             for(k=1;k<=Nz;k++)
                {
		     if(Lattice[i][j][k]==0)
		       {
                         if(count<=Nf)
                           {
                            wij[i][j][k]=gradw[count];
                            dens[i][j][k][0]=gradD[count];
                            if(strcmp(densread,"distinguished")==0)
                              {
                               dens1[i][j][k][0]=gradD1[count];
                               dens2[i][j][k][0]=gradD2[count];
                              }
                            sitesize[i][j][k]=gradp[count];
                            Ntot[i][j][k]=gradN[count];
                            pvol[i][j][k]=gradVpore[count];
                            vol_occ[i][j][k]=gradVocc[count];
                            inpore_c[i][j][k]=gradcc[count];
                            surface_pore[i][j][k]=gradpcc[count];
                            count+=2;
                           }
                         else
                           {
                            wij[i][j][k]=gradw[countt];
                            dens[i][j][k][0]=gradD[countt];
                            if(strcmp(densread,"distinguished")==0)
                              {
                               dens1[i][j][k][0]=gradD1[countt];
                               dens2[i][j][k][0]=gradD2[countt];
                              }
                            sitesize[i][j][k]=gradp[countt];
                            Ntot[i][j][k]=gradN[countt];
                            pvol[i][j][k]=gradVpore[countt];
                            vol_occ[i][j][k]=gradVocc[countt];
                            inpore_c[i][j][k]=gradcc[countt];
                            surface_pore[i][j][k]=gradpcc[countt];
                            countt-=2;
                           }
                       }  
                }
            }
         }
      }

     /* Gradiently organising property values in case of odd lattice dimensions */    
     if(Nx%2==1&&Ny%2==1&&Nz%2==1)
       {
        count=Nf;
        countt=0;
        for(i=((Nx-1)/2)+1;i>=1;i--)
           {
           for(j=1;j<=Ny;j++)
              {
               for(k=1;k<=Nz;k++)
                  {
                   if(Lattice[i][j][k]==0)
                      {
                       if(i==((Nx-1)/2)+1)
                            {
                             wij[i][j][k]=gradw[count];
                             dens[i][j][k][0]=gradD[count];
                             if(strcmp(densread,"distinguished")==0)
                               {
                                dens1[i][j][k][0]=gradD1[count];
                                dens2[i][j][k][0]=gradD2[count];
                               }
                             sitesize[i][j][k]=gradp[count];
                             Ntot[i][j][k]=gradN[count];
                             pvol[i][j][k]=gradVpore[count];
                             vol_occ[i][j][k]=gradVocc[count];
                             inpore_c[i][j][k]=gradcc[count];
                             surface_pore[i][j][k]=gradpcc[count];
                             count-=1;
                            }
                       if(i<((Nx-1)/2)+1)
                            {
                             wij[i][j][k]=gradw[count];
                             dens[i][j][k][0]=gradD[count];
                             if(strcmp(densread,"distinguished")==0)
                               {
                                dens1[i][j][k][0]=gradD1[count];
                                dens2[i][j][k][0]=gradD2[count];
                               }
                             sitesize[i][j][k]=gradp[count];
                             Ntot[i][j][k]=gradN[count];
                             pvol[i][j][k]=gradVpore[count];
                             vol_occ[i][j][k]=gradVocc[count];
                             inpore_c[i][j][k]=gradcc[count];
                             surface_pore[i][j][k]=gradpcc[count];
                             count-=1;
                             wij[i+2*countt][j][k]=gradw[count];
                             dens[i+2*countt][j][k][0]=gradD[count];
                             if(strcmp(densread,"distinguished")==0)
                               {
                                dens1[i+2*countt][j][k][0]=gradD1[count];
                                dens2[i+2*countt][j][k][0]=gradD2[count];
                               }
                             sitesize[i+2*countt][j][k]=gradp[count];
                             Ntot[i+2*countt][j][k]=gradN[count];
                             pvol[i+2*countt][j][k]=gradVpore[count];
                             vol_occ[i+2*countt][j][k]=gradVocc[count];
                             inpore_c[i+2*countt][j][k]=gradcc[count];
                             surface_pore[i+2*countt][j][k]=gradpcc[count];
                             count-=1;
                            }
                      }
                  }
              }
            countt++;
           }
       }     

    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                fprintf(outshifts,"%lf\n",wij[i][j][k]/(2.0*larmorfreq*PI));
               }
           }
       }
    for(k=1;k<=Nz;k++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(i=1;i<=Nx;i++)
               {
                fprintf(outgradpore,"%lf ",sitesize[i][j][k]);
               }
            fprintf(outgradpore,"\n");
           }
         fprintf(outgradpore,"---------------------------------------------\n");
        }
   fclose(outgradpore);
   fclose(outshifts);
  }

 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {    
	        for(l=0;l<=5;l++) 
		   {
		    if(alpha[i][j][k][l]==0.0) 
		      {
		       if(Enbar_mode==0)
		       {
		       alpha[i][j][k][l]=genergauss(meanEa,stdEa);
		       if(stdEa==0.0) alpha[i][j][k][l]=meanEa;
                       if(l==0)
                         {
                          i1=i+1;
                          if(i1>Nx) i1-=Nx;
                          alpha[i1][j][k][1]=alpha[i][j][k][l];
                         }
                       if(l==1)
                         {
                          i2=i-1;
                          if(i2<=0) i2+=Nx;
                          alpha[i2][j][k][0]=alpha[i][j][k][l];
                         }
                       if(l==2)
                         {
                          j1=j+1;
                          if(j1>Ny) j1-=Ny;
                          alpha[i][j1][k][3]=alpha[i][j][k][l];
                         }
                       if(l==3)
                         {
                          j2=j-1;
                          if(j2<=0) j2+=Ny;
                          alpha[i][j2][k][2]=alpha[i][j][k][l];
                         }
                       if(l==4)
                         {
                          k1=k+1;
                          if(k1>Nz) k1-=Nz;
                          alpha[i][j][k1][5]=alpha[i][j][k][l];
                         }
                       if(l==5)
                         {
                          k2=k-1;
                          if(k2<=0) k2+=Nz;
                          alpha[i][j][k2][4]=alpha[i][j][k][l];
                         }
                       }

                       if(Enbar_mode==1)
		       {
                       fitcst=0.0163494;
		       if(l==0)  
		         {
		          i1=i+1;
		          if(i1>Nx) i1-=Nx;
                          if(Lattice[i1][j][k]==0)
                             {  
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i1][j][k][1]=kB*T_ref*(inpore_c[i1][j][k]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i1][j][k]-log(fitcst)<0.0) alpha[i1][j][k][1]=0.0; 
                             }
                          if(Lattice[i1][j][k]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i1][j][k][1]=alpha[i][j][k][l];
                             }
		         }
		       if(l==1)  
		         {
		          i2=i-1;
		          if(i2<=0) i2+=Nx;
                          if(Lattice[i2][j][k]==0) 
                             { 
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i2][j][k][0]=kB*T_ref*(inpore_c[i2][j][k]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i2][j][k]-log(fitcst)<0.0) alpha[i2][j][k][0]=0.0;
                             } 
                          if(Lattice[i2][j][k]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i2][j][k][0]=alpha[i][j][k][l];
                             }
		         }
		       if(l==2)  
		         {
		          j1=j+1;
		          if(j1>Ny) j1-=Ny;
                          if(Lattice[i][j1][k]==0)
                             { 
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i][j1][k][3]=kB*T_ref*(inpore_c[i][j1][k]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i][j1][k]-log(fitcst)<0.0) alpha[i][j1][k][3]=0.0;
                             }  
                          if(Lattice[i][j1][k]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i][j1][k][3]=alpha[i][j][k][l];
                             }
		         }
		       if(l==3)  
		         {
		          j2=j-1;
		          if(j2<=0) j2+=Ny;
                          if(Lattice[i][j2][k]==0) 
                             { 
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i][j2][k][2]=kB*T_ref*(inpore_c[i][j2][k]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i][j2][k]-log(fitcst)<0.0) alpha[i][j2][k][2]=0.0;
                             } 
                          if(Lattice[i][j2][k]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i][j2][k][2]=alpha[i][j][k][l];
                             }
		         }
		       if(l==4)  
		         {
		          k1=k+1;
		          if(k1>Nz) k1-=Nz;
                          if(Lattice[i][j][k1]==0) 
                             { 
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i][j][k1][5]=kB*T_ref*(inpore_c[i][j][k1]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i][j][k1]-log(fitcst)<0.0) alpha[i][j][k1][5]=0.0;
                             } 
                          if(Lattice[i][j][k1]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i][j][k1][5]=alpha[i][j][k][l];
                             }
		         }
		       if(l==5)  
		         {
		          k2=k-1;
		          if(k2<=0) k2+=Nz;
                          if(Lattice[i][j][k2]==0) 
                             { 
                              alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                              alpha[i][j][k2][4]=kB*T_ref*(inpore_c[i][j][k2]-log(fitcst));
                              if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                              if(inpore_c[i][j][k2]-log(fitcst)<0.0) alpha[i][j][k2][4]=0.0;
                             } 
                          if(Lattice[i][j][k2]==1)
                             {
                              alpha[i][j][k][l]=pow(10,20);
                              alpha[i][j][k2][4]=alpha[i][j][k][l];
                             }
		         }
		      }
		  }
	       }
           }
       }
    }
 }

 for(i=1;i<=Nx;i++)
   {
    for(j=1;j<=Ny;j++)
       {
        for(k=1;k<=Nz;k++)
           {
            if(Lattice[i][j][k]==0)
              {
               for(l=0;l<=5;l++)
                  {
                   countb++;
                   fprintf(outbarr,"%d          %e\n",countb,alpha[i][j][k][l]);
                  }
              }
           }
       }
    }


 fclose(outxyz); 
 fclose(outpores); 
 fclose(outbarr); 
 fclose(outpsd);
 fclose(outpsdpart);

 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
		dens[i][j][k][0]/=part_func;
                Esite[i][j][k]=-kB*T*log(dens[i][j][k][0]);
                if(strcmp(densread,"distinguished")==0)
                  {
                   dens1[i][j][k][0]/=part_funct1;
                   E1site[i][j][k]=-kB*T*log(dens1[i][j][k][0]);
                   if(strcmp(solvent_activ,"solvent_ON")==0)
                     { 
                      dens2[i][j][k][0]/=part_funct2;
                      E2site[i][j][k]=-kB*T*log(dens2[i][j][k][0]);
                     }
                   if(strcmp(solvent_activ,"solvent_OFF")==0)
                     {
                      E2site[i][j][k]=0.0; 
                     }
                  }
	       }
	    }
	}
   }

 outenergy=fopen("poresize_energy.dat","w");
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
        {
         for(k=1;k<=Nz;k++)
            {
             if(Lattice[i][j][k]==0)
                {
                 fprintf(outenergy,"%lf		%lf\n",sitesize[i][j][k],Esite[i][j][k]);
                }
            }
        }
     }

  fclose(outenergy);

 out=fopen("Profiles.dat","w"); 
 fprintf(out,"# z - Energy - Frequency\n");
 for(k=1;k<=Nz;k++)
    {
     meanengz=0.0;
     meanfreqz=0.0;
     meandensz=0.0;
     zcount=0.0;
     for(j=1;j<=Ny;j++)
	{
	 for(i=1;i<=Nx;i++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
		meanengz+=Esite[i][j][k];
	        meanfreqz+=wij[i][j][k]/larmorfreq;
	  	meandensz+=dens[i][j][k][0];
		zcount++; 
	       } 
	    }
	}
     if(zcount!=0) {fprintf(out,"%e	%lf	%lf	%lf\n",k*a,meanengz/(zcount*1.0),meanfreqz/(zcount*1.0),meandensz/(zcount*1.0));}
   }
 
 fclose(out);
 
 printf("Reading OK!\n");
}


/* Equilibration function */
void equilibration()
{
 FILE *outdensvar;
 int i,j,k,l,m;
 int npore;
 int i1,i2,j1,j2,k1,k2;
 double parr,pleave;
 double parr1,parr2;
 double pleave1,pleave2;
 double ptr;
 char vardens[100];

 rho=fddmatrix(Nx+1,Ny+1,Nz+1,step+1);

 printf("Starting equilibration:\n"); 

 for(m=1;m<=step;m++)
    {
     for(i=1;i<=Nx;i++)
        {
         for(j=1;j<=Ny;j++)
            {
             for(k=1;k<=Nz;k++)
                {
                 if(Lattice[i][j][k]==0)
                   {
                    i1=i+1;
                    if(i1>Nx) i1-=Nx;
                    i2=i-1;
                    if(i2<=0) i2+=Nx;
                    j1=j+1;
                    if(j1>Ny) j1-=Ny;
                    j2=j-1;
                    if(j2<=0) j2+=Ny;
                    k1=k+1;
                    if(k1>Nz) k1-=Nz;
                    k2=k-1;
                    if(k2<=0) k2+=Nz;

                    /* Prob of arriving at site ijk */
	    	    parr=0.0;
                    parr1=0.0;
                    parr2=0.0;
		    ptr=exp(-alpha[i1][j][k][1]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i1][j][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i1][j][k][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i1][j][k])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i1][j][k][0];
                         }
                       if(E2site[i][j][k]<=E2site[i1][j][k])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i1][j][k][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i1][j][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i1][j][k])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i1][j][k][0]*exp(-(E1site[i][j][k]-E1site[i1][j][k])/(kB*T))*ptr;
                         }
                       if(E2site[i][j][k]>E2site[i1][j][k])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i1][j][k][0]*exp(-(E2site[i][j][k]-E2site[i1][j][k])/(kB*T))*ptr;
                         }
                      }
		    ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i2][j][k])
                      { 
                       parr+=(1/(dir*1.0))*dens[i2][j][k][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i2][j][k])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i2][j][k][0];
                         }
                       if(E2site[i][j][k]<=E2site[i2][j][k])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i2][j][k][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i2][j][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i2][j][k])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i2][j][k][0]*exp(-(E1site[i][j][k]-E1site[i2][j][k])/(kB*T))*ptr;
                         }
                       if(E2site[i][j][k]>E2site[i2][j][k])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i2][j][k][0]*exp(-(E2site[i][j][k]-E2site[i2][j][k])/(kB*T))*ptr;
                         }
                      }
		    ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i][j1][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i][j1][k][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i][j1][k])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i][j1][k][0];
                         }
                       if(E2site[i][j][k]<=E2site[i][j1][k])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i][j1][k][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i][j1][k])
                      { 
                       parr+=(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i][j1][k])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i][j1][k][0]*exp(-(E1site[i][j][k]-E1site[i][j1][k])/(kB*T))*ptr;
                         }
                       if(E2site[i][j][k]>E2site[i][j1][k])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i][j1][k][0]*exp(-(E2site[i][j][k]-E2site[i][j1][k])/(kB*T))*ptr;
                         }
                      }
		    ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i][j2][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i][j2][k][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i][j2][k])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i][j2][k][0];
                         }
                       if(E2site[i][j][k]<=E2site[i][j2][k])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i][j2][k][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i][j2][k]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i][j2][k])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i][j2][k][0]*exp(-(E1site[i][j][k]-E1site[i][j2][k])/(kB*T))*ptr;
                         }
                       if(E2site[i][j][k]>E2site[i][j2][k])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i][j2][k][0]*exp(-(E2site[i][j][k]-E2site[i][j2][k])/(kB*T))*ptr;
                         }
                      }
		    ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i][j][k1]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i][j][k1][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i][j][k1])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k1][0]; 
                         }
                       if(E2site[i][j][k]<=E2site[i][j][k1])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k1][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i][j][k1])
                      { 
                       parr+=(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i][j][k1])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i][j][k1][0]*exp(-(E1site[i][j][k]-E1site[i][j][k1])/(kB*T))*ptr;
                         }
                       if(E2site[i][j][k]>E2site[i][j][k1])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i][j][k1][0]*exp(-(E2site[i][j][k]-E2site[i][j][k1])/(kB*T))*ptr;
                         }
                      }
		    ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	    	    if(Esite[i][j][k]<=Esite[i][j][k2])
                      { 
                       parr+=(1/(dir*1.0))*dens[i][j][k2][0]*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]<=E1site[i][j][k2])
                         {
                          parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k2][0]; 
                         }
                       if(E2site[i][j][k]<=E2site[i][j][k2])
                         {
                          parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k2][0];
                         }
                      }
	    	    if(Esite[i][j][k]>Esite[i][j][k2]) 
                      {
                       parr+=(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
                      }
                    if(strcmp(densread,"distinguished")==0)
                      {
                       if(E1site[i][j][k]>E1site[i][j][k2])
                         {
                          parr1+=(1/(dir*1.0))*dens1[i][j][k2][0]*exp(-(E1site[i][j][k]-E1site[i][j][k2])/(kB*T))*ptr;
                         }   
                       if(E2site[i][j][k]>E2site[i][j][k2])
                         {
                          parr2+=(1/(dir*1.0))*dens2[i][j][k2][0]*exp(-(E2site[i][j][k]-E2site[i][j][k2])/(kB*T))*ptr;
                         }
                      }

                    /* Prob of leaving the site ijk */
	    	    pleave=0.0;
                    pleave1=0.0;
                    pleave2=0.0;
	    	    if(Lattice[i1][j][k]==0)
	      	      {
                       ptr=exp(-alpha[i][j][k][0]/(kB*T));
		       if(Esite[i1][j][k]<=Esite[i][j][k])
                         { 
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i1][j][k]<=E1site[i][j][k])
                            {
                             pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                            }
                          if(E2site[i1][j][k]<=E2site[i][j][k])
                            {
                             pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                            }
                         }
		       if(Esite[i1][j][k]>Esite[i][j][k]) 
                         {
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i1][j][k]>E1site[i][j][k])
                            {
                             pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i1][j][k]-E1site[i][j][k])/(kB*T))*ptr;
                            }
                          if(E2site[i1][j][k]>E2site[i][j][k])
                            {
                             pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i1][j][k]-E2site[i][j][k])/(kB*T))*ptr;
                            }
                         }
	     	      }   
	    	    if(Lattice[i2][j][k]==0)
	      	      {
                       ptr=exp(-alpha[i][j][k][1]/(kB*T));
		       if(Esite[i2][j][k]<=Esite[i][j][k])
                         { 
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i2][j][k]<=E1site[i][j][k])
                            {
                             pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                            }
                          if(E2site[i2][j][k]<=E2site[i][j][k])
                            {
                             pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                            }
                         }
		       if(Esite[i2][j][k]>Esite[i][j][k])
                         { 
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i2][j][k]>E1site[i][j][k])
                            {
                             pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i2][j][k]-E1site[i][j][k])/(kB*T))*ptr;
                            }
                          if(E2site[i2][j][k]>E2site[i][j][k])
                            {
                             pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i2][j][k]-E2site[i][j][k])/(kB*T))*ptr;
                            }
                         }
	      	      }
	    	    if(Lattice[i][j1][k]==0)
	      	      {
                       ptr=exp(-alpha[i][j][k][2]/(kB*T));
		       if(Esite[i][j1][k]<=Esite[i][j][k]) 
                         {
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i][j1][k]<=E1site[i][j][k])
                            {  
                             pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                            }
                          if(E2site[i][j1][k]<=E2site[i][j][k])
                            {
                             pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                            }
                         }
		       if(Esite[i][j1][k]>Esite[i][j][k])
                         { 
                          pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr;
                         }
                       if(strcmp(densread,"distinguished")==0)
                         {
                          if(E1site[i][j1][k]>E1site[i][j][k])
                            {
                             pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j1][k]-E1site[i][j][k])/(kB*T))*ptr;
                            }
                          if(E2site[i][j1][k]>E2site[i][j][k])
                            {
                             pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j1][k]-E2site[i][j][k])/(kB*T))*ptr;
                            }
                         }
	      	      }
	    	   if(Lattice[i][j2][k]==0)
	      	     {
                      ptr=exp(-alpha[i][j][k][3]/(kB*T));
		      if(Esite[i][j2][k]<=Esite[i][j][k])
                        { 
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j2][k]<=E1site[i][j][k])
                           {
                            pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                           }
                         if(E2site[i][j2][k]<=E2site[i][j][k])
                           {
                            pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                           }
                        }
		      if(Esite[i][j2][k]>Esite[i][j][k]) 
                        {
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j2][k]>E1site[i][j][k])
                           {
                            pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j2][k]-E1site[i][j][k])/(kB*T))*ptr;
                           }
                         if(E2site[i][j2][k]>E2site[i][j][k])
                           {
                            pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j2][k]-E2site[i][j][k])/(kB*T))*ptr;
                           }
                        }
	      	     }
	    	   if(Lattice[i][j][k1]==0)
	      	     {
                      ptr=exp(-alpha[i][j][k][4]/(kB*T));
		      if(Esite[i][j][k1]<=Esite[i][j][k]) 
                        {
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j][k1]<=E1site[i][j][k])
                           {
                            pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                           }
                         if(E2site[i][j][k1]<=E2site[i][j][k])
                           {
                            pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                           }
                        }
		      if(Esite[i][j][k1]>Esite[i][j][k]) 
                        {
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j][k1]>E1site[i][j][k])
                           {
                            pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k1]-E1site[i][j][k])/(kB*T))*ptr;
                           }
                         if(E2site[i][j][k1]>E2site[i][j][k])
                           {
                            pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k1]-E2site[i][j][k])/(kB*T))*ptr;
                           }
                        }
	      	     }
	    	   if(Lattice[i][j][k2]==0)
	      	     {
                      ptr=exp(-alpha[i][j][k][5]/(kB*T));
		      if(Esite[i][j][k2]<=Esite[i][j][k]) 
                        {
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j][k2]<=E1site[i][j][k])
                           {
                            pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                           }
                         if(E2site[i][j][k2]<=E2site[i][j][k])
                           {
                            pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                           }
                        }
		      if(Esite[i][j][k2]>Esite[i][j][k])
                        { 
                         pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr;
                        }
                      if(strcmp(densread,"distinguished")==0)
                        {
                         if(E1site[i][j][k2]>E1site[i][j][k])
                           {
                            pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k2]-E1site[i][j][k])/(kB*T))*ptr;
                           }
                         if(E2site[i][j][k2]>E2site[i][j][k])
                           {
                            pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k2]-E2site[i][j][k])/(kB*T))*ptr;
                           }
                        }
	      	     }

                   /* Evolution of the density */
                   dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
                   rho[i][j][k][m]=dens[i][j][k][1];
                   if(strcmp(densread,"distinguished")==0)
                     {
                      dens1[i][j][k][1]=dens1[i][j][k][0]+parr1-pleave1;
                      dens2[i][j][k][1]=dens2[i][j][k][0]+parr2-pleave2;
                      if(strcmp(iontype,"anion")==0)
                        {
                         vol_occ[i][j][k]=Vanion*dens[i][j][k][1]*part_funct+Vcation*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                         inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                        }
                      if(strcmp(iontype,"cation")==0)
                        {
                         vol_occ[i][j][k]=Vcation*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                         inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                        }
                      if(strcmp(iontype,"solvent")==0)
                        {
                         vol_occ[i][j][k]=Vsolvent*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vcation*dens2[i][j][k][1]*part_funct2;
                         inpore_c[i][j][k]=(((dens1[i][j][k][1]*part_funct1+dens2[i][j][k][1]*part_funct2)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                        }
                     }
                   if(strcmp(densread,"undistinguished")==0)
                     {
                      vol_occ[i][j][k]=Vspec*dens[i][j][k][1]*part_funct;
                      inpore_c[i][j][k]=((dens[i][j][k][1]*part_funct*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }

                   dens[i][j][k][0]=dens[i][j][k][1];
                   if(strcmp(densread,"distinguished")==0)
                     {
                      dens1[i][j][k][0]=dens1[i][j][k][1];
                      dens2[i][j][k][0]=dens2[i][j][k][1];
                     }
                  }
               }
            }
        }

      if(Enbar_mode==1)
       {	      
       for(i=1;i<=Nx;i++)
          {
           for(j=1;j<=Ny;j++)
              {
               for(k=1;k<=Nz;k++)
                  {
                   if(Lattice[i][j][k]==0)
                     {
                      for(l=0;l<=5;l++)
                         {
                          if(alpha[i][j][k][l]!=pow(10,20))
                            {
                             alpha[i][j][k][l]=kB*T_ref*(inpore_c[i][j][k]-log(fitcst));
                             if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                            }
                         }
                     }
                  }
              }
          }
       }
    }

   if(strcmp(wrtdns,"wrtdns_on")==0)
     {
      system("rm -rf Pore-density-profile");
      system("mkdir Pore-density-profile");
      for(i=1;i<=Nx;i++)
         {
          for(j=1;j<=Ny;j++)
             {
              for(k=1;k<=Nz;k++)
                 {
                  if(Lattice[i][j][k]==0)
                    {
                     sprintf(vardens,"./Pore-density-profile/Density_variation_pore_%d%d%d.dat",i,j,k);
                     outdensvar=fopen(vardens,"w");
                     for(m=1;m<=step;m=m+nsmple)
                        {
                         fprintf(outdensvar,"%lf     %lf\n",m*dtime,rho[i][j][k][m]);
                        } 
                        fclose(outdensvar);
                    }
                 }
             }
        }
     }                    

    /* Verifying equilibrium of pore densities */
    npore=0.0; 
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  { 
                   if(fabs(rho[i][j][k][step]-rho[i][j][k][step-1])<thrsld)
                     {
                      if(fabs(rho[i][j][k][step]-rho[i][j][k][step-2])<thrsld)
                        {
                         if(fabs(rho[i][j][k][step]-rho[i][j][k][step-3])<thrsld)  
                           {
                            npore++;
                           }
                        }
                     }
                  }
               }
           }
       }
    if(npore==Nf)
      {  
       printf("Equilibration process succeeded!");
      }
    else
      {
       printf("%lf %% of sites are not in a stationary state!",(npore*100.0)/Nf);
      }
}


/* Main function: propagation of the msd and NMR signal */
void propagation()
{
 FILE *out,*outx,*outy,*outz;
 FILE *outdenstot,*outG,*outGsmp;
 int i,j,k,l,t;
 int i1,i2,j1,j2,k1,k2,treal;
 double parr,pleave,denstot;
 double ptr,locpfunc;
 double vxmoy,vymoy,vzmoy,Dlargest,Dsmallest,anisotropy;
 double magnitude,phase; 
 double newRe,newIm;

 FILE *outDt,*outmsd;
 double msd,msdx,msdy,msdz;
 double summsd1,summsd2;
 double summsd1x,summsd2x;
 double summsd1y,summsd2y;
 double summsd1z,summsd2z;
 double msdi,msdxi,msdyi,msdzi;
 double msdf,msdxf,msdyf,msdzf;
 double Df,Dfx,Dfy,Dfz;
 double parr1,parr2,pleave1,pleave2;

 /* locpfunc is the "local partition function" */

 outdenstot=fopen("density_tot.out","w");
 outG=fopen("Signal.dat","w");
 fprintf(outG,"# Time - Real part - Imaginary part - Magnitude - Phase\n");
 outGsmp=fopen("Signal_smp.dat","w");
 fprintf(outGsmp,"# Time - Real part - Imaginary part - Magnitude - Phase\n");

 out=fopen("vacf.out","w");
 outx=fopen("vacfx.out","w");
 outy=fopen("vacfy.out","w");
 outz=fopen("vacfz.out","w");
 outmsd=fopen("MSD.dat","w");
 outDt=fopen("time_dep_Diffusion_coeff-D_over_D0.dat","w");
 
 /* 1st term of the velocity autocorrelation function = mean of v(0)^2/2. */
 denstot=0.0;
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if((Lattice[i][j][k]==0)&&(dens[i][j][k][0]!=0))
	       {
	       /* Neighbouring sites and application of PBC */
		i1=i+1;
		if(i1>Nx) i1-=Nx;
		i2=i-1;	
		if(i2<=0) i2+=Nx;
		j1=j+1;
		if(j1>Ny) j1-=Ny;
		j2=j-1;
		if(j2<=0) j2+=Ny;
		k1=k+1;
		if(k1>Nz) k1-=Nz;
		k2=k-1;
		if(k2<=0) k2+=Nz;
		/* ptr is the probability of the transition */
		if(Lattice[i1][j][k]==0) 
		  {
                   ptr=1.0*exp(-alpha[i][j][k][0]/(kB*T));
                   if(Esite[i1][j][k]>Esite[i][j][k]) ptr*=exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T));
		   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfx+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfx[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfx[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfx[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 
		if(Lattice[i2][j][k]==0) 
		  {
                   ptr=1.0*exp(-alpha[i][j][k][1]/(kB*T));
                   if(Esite[i2][j][k]>Esite[i][j][k]) ptr*=exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T));
		   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfx+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfx[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfx[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfx[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 	
		if(Lattice[i][j1][k]==0) 
		  {
                   ptr=1.0*exp(-alpha[i][j][k][2]/(kB*T));
                   if(Esite[i][j1][k]>Esite[i][j][k]) ptr*=exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T));
		   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfy+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfy[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfy[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfy[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 
		if(Lattice[i][j2][k]==0) 
		  {
                   ptr=1.0*exp(-alpha[i][j][k][3]/(kB*T));
                   if(Esite[i][j2][k]>Esite[i][j][k]) ptr*=exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T));
		   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfy+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfy[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfy[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfy[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 
		if(Lattice[i][j][k1]==0) 	
		  {
                   ptr=1.0*exp(-alpha[i][j][k][4]/(kB*T));
                   if(Esite[i][j][k1]>Esite[i][j][k]) ptr*=exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T));
	 	   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfz+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfz[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfz[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfz[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 
		if(Lattice[i][j][k2]==0) 
		  {
                   ptr=1.0*exp(-alpha[i][j][k][5]/(kB*T));
                   if(Esite[i][j][k2]>Esite[i][j][k]) ptr*=exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T));
		   lattvacf+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; lattvacfz+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; 
		   xvacf[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; xvacfz[i]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   yvacf[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; yvacfz[j]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		   zvacf[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr; zvacfz[k]+=(1/(dir*1.0))*dens[i][j][k][0]*v0*v0*ptr;
		  } 
	   }
	 denstot+=dens[i][j][k][0];
	}
    }
 }

 fprintf(outdenstot,"%e	%lf\n",0.0*dtime,denstot);

 /* For the calculation of diffusion coefficients */
 fprintf(out,"%e %e\n",0.0,lattvacf);
 fprintf(outx,"%e %e\n",0.0,lattvacfx);
 fprintf(outy,"%e %e\n",0.0,lattvacfy);
 fprintf(outz,"%e %e\n",0.0,lattvacfz);
 sumvacf+=lattvacf*dtime/(2.0*denstot);
 sumvacfx+=lattvacfx*dtime/(2.0*denstot);
 sumvacfy+=lattvacfy*dtime/(2.0*denstot);
 sumvacfz+=lattvacfz*dtime/(2.0*denstot);
 for(i=1;i<=Nx;i++)
    {xsumvacf[i]+=xvacf[i]*dtime/2.0; xsumvacfx[i]+=xvacfx[i]*dtime/2.0; xsumvacfy[i]+=xvacfy[i]*dtime/2.0; xsumvacfz[i]+=xvacfz[i]*dtime/2.0;}
 for(j=1;j<=Ny;j++)
    {ysumvacf[j]+=yvacf[j]*dtime/2.0; ysumvacfx[j]+=yvacfx[j]*dtime/2.0; ysumvacfy[j]+=yvacfy[j]*dtime/2.0; ysumvacfz[j]+=yvacfz[j]*dtime/2.0;}
 for(k=1;k<=Nz;k++)
    {zsumvacf[k]+=zvacf[k]*dtime/2.0; zsumvacfx[k]+=zvacfx[k]*dtime/2.0; zsumvacfy[k]+=zvacfy[k]*dtime/2.0; zsumvacfz[k]+=zvacfz[k]*dtime/2.0;}


 msdvacf[0]=lattvacf*dtime/denstot;
 msdvacfx[0]=lattvacfx*dtime/denstot;
 msdvacfy[0]=lattvacfy*dtime/denstot;
 msdvacfz[0]=lattvacfz*dtime/denstot;
 fprintf(outDt,"%e      %lf\n",0.0*dtime,(sumvacf/(dim*1.0))/(a*a/(2.0*dim*dtime)));

 
 /* At t=0, G(0)=1 if no spin echo before */
 ReGtot[0]=1.0;
 ImGtot[0]=0.0;
 ReGtotsmp[0]=1.0;
 ImGtotsmp[0]=0.0;
 magnitude=sqrt(ReGtot[0]*ReGtot[0]+ImGtot[0]*ImGtot[0]);
 phase=atan(ImGtot[0]/ReGtot[0]);
 fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",0.0,ReGtot[0],ImGtot[0],magnitude,phase);
 fprintf(outGsmp,"%e	%lf 	%lf	%lf	%lf\n",0.0,ReGtotsmp[0],ImGtotsmp[0],magnitude,phase);

 /* Initiation of the propagated function which is
 the probability of arriving at site (i,j) at time 1 
 weighted by the initial velocity v0 */
 denstot=0.0;
 lattvacf=0.0;
 lattvacfx=0.0;
 lattvacfy=0.0;
 lattvacfz=0.0;

 for(i=1;i<=Nx;i++)
    {xvacf[i]=0.0; xvacfx[i]=0.0; xvacfy[i]=0.0; xvacfz[i]=0.0;}
 for(j=1;j<=Ny;j++)
    {yvacf[j]=0.0; yvacfx[j]=0.0; yvacfy[j]=0.0; yvacfz[j]=0.0;}
 for(k=1;k<=Nz;k++)
    {zvacf[k]=0.0; zvacfx[k]=0.0; zvacfy[k]=0.0; zvacfz[k]=0.0;}
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     vxmoy=0.0; vymoy=0.0; vzmoy=0.0;
	     if(Lattice[i][j][k]==0)
	       {
		/* Neighbouring sites and application of PBC */
		i1=i+1;
		if(i1>Nx) i1-=Nx;
		i2=i-1;
		if(i2<=0) i2+=Nx;
		j1=j+1;
		if(j1>Ny) j1-=Ny;
		j2=j-1;
		if(j2<=0) j2+=Ny;
		k1=k+1;
		if(k1>Nz) k1-=Nz;
		k2=k-1;
		if(k2<=0) k2+=Nz;

                pVtot[i][j][k]=vol_occ[i][j][k];
                pinpore_c[i][j][k]=inpore_c[i][j][k];

	    	/* Average velocity in site (i,j) */
	    	if((Lattice[i1][j][k]==0)&&(Esite[i1][j][k]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   vxmoy=ptr*(1/(dir*1.0))*v0;
	  	  }
	    	if((Lattice[i1][j][k]==0)&&(Esite[i1][j][k]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   vxmoy=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i1][j][k]-Esite[i][j][k]));
		  }
	    	if((Lattice[i2][j][k]==0)&&(Esite[i2][j][k]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
		   vxmoy-=ptr*(1/(dir*1.0))*v0;
		  }
	    	if((Lattice[i2][j][k]==0)&&(Esite[i2][j][k]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
		   vxmoy-=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i2][j][k]-Esite[i][j][k]));
		  }
	    	if((Lattice[i][j1][k]==0)&&(Esite[i][j1][k]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
		   vymoy=ptr*(1/(dir*1.0))*v0;
		  }
	    	if((Lattice[i][j1][k]==0)&&(Esite[i][j1][k]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
		   vymoy=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i][j1][k]-Esite[i][j][k]));
		  }
	    	if((Lattice[i][j2][k]==0)&&(Esite[i][j2][k]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
		   vymoy-=ptr*(1/(dir*1.0))*v0;
		  }
	    	if((Lattice[i][j2][k]==0)&&(Esite[i][j2][k]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
		   vymoy-=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i][j2][k]-Esite[i][j][k]));
		  }
	    	if((Lattice[i][j][k1]==0)&&(Esite[i][j][k1]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
		   vzmoy=ptr*(1/(dir*1.0))*v0;
		  }
	    	if((Lattice[i][j][k1]==0)&&(Esite[i][j][k1]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
		   vzmoy=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i][j][k1]-Esite[i][j][k]));
		  }
	    	if((Lattice[i][j][k2]==0)&&(Esite[i][j][k2]<=Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
		   vzmoy-=ptr*(1/(dir*1.0))*v0;
		  }
	    	if((Lattice[i][j][k2]==0)&&(Esite[i][j][k2]>Esite[i][j][k])) 
		  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
		   vzmoy-=ptr*(1/(dir*1.0))*v0*exp(-(Esite[i][j][k2]-Esite[i][j][k]));
		  }
	    	/* Probability of arriving in site i,j at t=1 weighted by the velocity */
	    	locpfunc=0.0; 
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i1][j][k][1]/(kB*T));
                   /* what arrives */
	           probv0[i][j][k][0][1]-=dens[i1][j][k][0]*v0*ptr; 
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* what stays */
		   /* if it stays, velocity equal zero */
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	          }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i1][j][k][1]/(kB*T));
                   /* what arrives */
	           probv0[i][j][k][0][1]-=dens[i1][j][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* what stays */
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	       	   probv0[i][j][k][0][1]+=dens[i2][j][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i2][j][k][0]*ptr;
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	       	   probv0[i][j][k][0][1]+=dens[i2][j][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	       	   probv0[i][j][k][1][1]-=dens[i][j1][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j1][k][0]*ptr;
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	       	   probv0[i][j][k][1][1]-=dens[i][j1][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	       	   probv0[i][j][k][1][1]+=dens[i][j2][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j2][k][0]*ptr;
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	 }
	    	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	       	   probv0[i][j][k][1][1]+=dens[i][j2][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  }
	    	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	       	   probv0[i][j][k][2][1]-=dens[i][j][k1][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j][k1][0]*ptr;
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	       	   probv0[i][j][k][2][1]-=dens[i][j][k1][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	       	   probv0[i][j][k][2][1]+=dens[i][j][k2][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j][k2][0]*ptr;
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);

	      	 }
	    	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	       	   probv0[i][j][k][2][1]+=dens[i][j][k2][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  }
	     	locpfunc+=(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2])*dens[i][j][k][0];
	        for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t=1 */
	    	parr=0.0;
                parr1=0.0;
                parr2=0.0;
                ptr=exp(-alpha[i1][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i1][j][k][0];
                  }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j][k]<=E1site[i1][j][k])
                       {
                        parr1+=ptr*(1/(dir*1.0))*dens1[i1][j][k][0];
                       }
                     if(E2site[i][j][k]<=E2site[i1][j][k])
                       {
                        parr2+=ptr*(1/(dir*1.0))*dens2[i1][j][k][0];
                       }
                    }
	    	if(Esite[i][j][k]>Esite[i1][j][k]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i1][j][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i1][j][k][0]*exp(-(E1site[i][j][k]-E1site[i1][j][k])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i1][j][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i1][j][k][0]*exp(-(E2site[i][j][k]-E2site[i1][j][k])/(kB*T));
                     }
                  }
                ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i2][j][k][0];
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i2][j][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i2][j][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i2][j][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i2][j][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i2][j][k]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i2][j][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i2][j][k][0]*exp(-(E1site[i][j][k]-E1site[i2][j][k])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i2][j][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i2][j][k][0]*exp(-(E2site[i][j][k]-E2site[i2][j][k])/(kB*T));
                     }
                  }
		ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i][j1][k][0];
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j1][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j1][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j1][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j1][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j1][k]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j1][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j1][k][0]*exp(-(E1site[i][j][k]-E1site[i][j1][k])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i][j1][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j1][k][0]*exp(-(E2site[i][j][k]-E2site[i][j1][k])/(kB*T));
                     }
                  }
		ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i][j2][k][0];
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j2][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j2][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j2][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j2][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j2][k])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j2][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j2][k][0]*exp(-(E1site[i][j][k]-E1site[i][j2][k])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i][j2][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j2][k][0]*exp(-(E2site[i][j][k]-E2site[i][j2][k])/(kB*T));
                     }
                  }
		ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i][j][k1][0];
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j][k1])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k1][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j][k1])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k1][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j][k1])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j][k1])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k1][0]*exp(-(E1site[i][j][k]-E1site[i][j][k1])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i][j][k1])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k1][0]*exp(-(E2site[i][j][k]-E2site[i][j][k1])/(kB*T));
                     }
                  }
		ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2])
                  { 
                   parr+=ptr*(1/(dir*1.0))*dens[i][j][k2][0];
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j][k2])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k2][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j][k2])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k2][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j][k2]) 
                  {
                   parr+=ptr*(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T));
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j][k2])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k2][0]*exp(-(E1site[i][j][k]-E1site[i][j][k2])/(kB*T));
                     }
                   if(E2site[i][j][k]>E2site[i][j][k2])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k2][0]*exp(-(E2site[i][j][k]-E2site[i][j][k2])/(kB*T));
                     }
                  }
	    	/* Probability of leaving site i,j at t=1 */
	    	pleave=0.0;
                pleave1=0.0;
                pleave2=0.0;
	    	if(Lattice[i1][j][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   if(Esite[i1][j][k]<=Esite[i][j][k])
                     { 
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i1][j][k]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i1][j][k]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i1][j][k]>Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i1][j][k]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i1][j][k]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i1][j][k]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i1][j][k]-E2site[i][j][k])/(kB*T));
                        }
                     }
	      	  }
	    	if(Lattice[i2][j][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][1]/(kB*T));
		   if(Esite[i2][j][k]<=Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i2][j][k]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i2][j][k]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i2][j][k]>Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i2][j][k]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i2][j][k]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i2][j][k]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i2][j][k]-E2site[i][j][k])/(kB*T));
                        }
                     }
		  }
	    	if(Lattice[i][j1][k]==0)
	          {
		   ptr=exp(-alpha[i][j][k][2]/(kB*T));
		   if(Esite[i][j1][k]<=Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j1][k]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i][j1][k]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i][j1][k]>Esite[i][j][k])
                     { 
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j1][k]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j1][k]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i][j1][k]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j1][k]-E2site[i][j][k])/(kB*T));
                        }
                     }
	          }
	    	if(Lattice[i][j2][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][3]/(kB*T));
		   if(Esite[i][j2][k]<=Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j2][k]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i][j2][k]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i][j2][k]>Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j2][k]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j2][k]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i][j2][k]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j2][k]-E2site[i][j][k])/(kB*T));
                        }
                     }
	          }
	    	if(Lattice[i][j][k1]==0)
	          {
		   ptr=exp(-alpha[i][j][k][4]/(kB*T));
		   if(Esite[i][j][k1]<=Esite[i][j][k])
                     { 
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j][k1]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i][j][k1]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i][j][k1]>Esite[i][j][k])
                     { 
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j][k1]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k1]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i][j][k1]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k1]-E2site[i][j][k])/(kB*T));
                        }
                     }
	          }
	    	if(Lattice[i][j][k2]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][5]/(kB*T));
		   if(Esite[i][j][k2]<=Esite[i][j][k])
                     { 
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j][k2]<=E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                        }
                      if(E2site[i][j][k2]<=E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                        }
                     }
		   if(Esite[i][j][k2]>Esite[i][j][k]) 
                     {
                      pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T));
                     }
                   if(strcmp(densread,"distinguished")==0)
                     {
                      if(E1site[i][j][k2]>E1site[i][j][k])
                        {
                         pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k2]-E1site[i][j][k])/(kB*T));
                        }
                      if(E2site[i][j][k2]>E2site[i][j][k])
                        {
                         pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k2]-E2site[i][j][k])/(kB*T));
                        }
                     }
	          }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
                if(strcmp(densread,"distinguished")==0)
                  {
                   dens1[i][j][k][1]=dens1[i][j][k][0]+parr1-pleave1;
                   dens2[i][j][k][1]=dens2[i][j][k][0]+parr2-pleave2;
                   if(strcmp(iontype,"anion")==0)
                     {
                      vol_occ[i][j][k]=Vanion*dens[i][j][k][1]*part_funct+Vcation*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                   if(strcmp(iontype,"cation")==0)
                     {
                      vol_occ[i][j][k]=Vcation*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                   if(strcmp(iontype,"solvent")==0)
                     {
                      vol_occ[i][j][k]=Vsolvent*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vcation*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens1[i][j][k][1]*part_funct1+dens2[i][j][k][1]*part_funct2)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                  }
                if(strcmp(densread,"undistinguished")==0)
                  {
                   vol_occ[i][j][k]=Vspec*dens[i][j][k][1]*part_funct;
                   inpore_c[i][j][k]=((dens[i][j][k][1]*part_funct*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                  }
	    	/* Initiation of the function (complex function) allowing for the calculation of the spectra */
	    	ReG[i][j][k][1]=cos(wij[i][j][k]*dtime);		/* Real part */
	    	ImG[i][j][k][1]=sin(wij[i][j][k]*dtime);		/* Imaginary part */
	    	ReGtot[1]+=dens[i][j][k][1]*ReG[i][j][k][1];
	    	ImGtot[1]+=dens[i][j][k][1]*ImG[i][j][k][1];
	       }
	     lattvacf+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; lattvacf+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     lattvacf+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     lattvacfx+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; lattvacfy+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     lattvacfz+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     xvacfx[i]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; xvacfy[i]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     xvacfz[i]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     yvacfx[j]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; yvacfy[j]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     yvacfz[j]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     zvacfx[k]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; zvacfy[k]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	     zvacfz[k]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	     denstot+=dens[i][j][k][1];
	    }
	}
    }

 if(Enbar_mode==1)
   {
    for(i=1;i<=Nx;i++)
       {
        for(j=1;j<=Ny;j++)
           {
            for(k=1;k<=Nz;k++)
               {
                if(Lattice[i][j][k]==0)
                  {
                   for(l=0;l<=5;l++)
                      {
                       if(alpha[i][j][k][l]!=pow(10,20))
                         {
                          alpha[i][j][k][l]=alpha[i][j][k][l]/(kB*T_ref*(pinpore_c[i][j][k]-log(fitcst)));
                          alpha[i][j][k][l]=alpha[i][j][k][l]*(kB*T_ref*(inpore_c[i][j][k]-log(fitcst)));
                          if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                         }
                      }
                  }
               }
           }
       }
   } 
 
 fprintf(outdenstot,"%e	%lf\n",1.0*dtime,denstot);

 ReGtot[1]/=denstot;
 ImGtot[1]/=denstot;
 magnitude=sqrt(ReGtot[1]*ReGtot[1]+ImGtot[1]*ImGtot[1]);
 phase=atan(ImGtot[1]/ReGtot[1]);
 fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",dtime,ReGtot[1],ImGtot[1],magnitude,phase);
 if(nsample==1)
   {
    ReGtotsmp[1]=ReGtot[1];
    ImGtotsmp[1]=ImGtot[1];
    fprintf(outGsmp,"%e	%lf 	%lf	%lf	%lf\n",dtime,ReGtotsmp[1],ImGtotsmp[1],magnitude,phase);
    treal=2;
   }
 if(nsample>1) treal=1;
 fprintf(out,"%e %e\n",dtime,lattvacf);
 fprintf(outx,"%e %e\n",dtime,lattvacfx);
 fprintf(outy,"%e %e\n",dtime,lattvacfy);
 fprintf(outz,"%e %e\n",dtime,lattvacfz);
 sumvacf+=lattvacf*dtime/denstot;
 sumvacfx+=lattvacfx*dtime/denstot;
 sumvacfy+=lattvacfy*dtime/denstot;
 sumvacfz+=lattvacfz*dtime/denstot;
 for(i=1;i<=Nx;i++)
    {xsumvacf[i]+=xvacf[i]*dtime; xsumvacfx[i]+=xvacfx[i]*dtime; xsumvacfy[i]+=xvacfy[i]*dtime; xsumvacfz[i]+=xvacfz[i]*dtime;}
 for(j=1;j<=Ny;j++)
    {ysumvacf[j]+=yvacf[j]*dtime; ysumvacfx[j]+=yvacfx[j]*dtime; ysumvacfy[j]+=yvacfy[j]*dtime; ysumvacfz[j]+=yvacfz[j]*dtime;}
 for(k=1;k<=Nz;k++)
    {zsumvacf[k]+=zvacf[k]*dtime; zsumvacfx[k]+=zvacfx[k]*dtime; zsumvacfy[k]+=zvacfy[k]*dtime; zsumvacfz[k]+=zvacfz[k]*dtime;}

 msdvacf[1]=lattvacf*dtime/denstot;
 msdvacfx[1]=lattvacfx*dtime/denstot;
 msdvacfy[1]=lattvacfy*dtime/denstot;
 msdvacfz[1]=lattvacfz*dtime/denstot;
 fprintf(outDt,"%e      %lf\n",1.0*dtime,(sumvacf/(dim*1.0))/(a*a/(2.0*dim*dtime)));

 printf("\n");

 /* Propagation */

 for(t=2;t<=Nmin;t++)
 {
  lattvacf=0.0;
  lattvacfx=0.0;
  lattvacfy=0.0;
  lattvacfz=0.0;


  for(i=1;i<=Nx;i++) {xvacf[i]=0.0; xvacfx[i]=0.0; xvacfy[i]=0.0; xvacfz[i]=0.0;}
  for(j=1;j<=Ny;j++) {yvacf[j]=0.0; yvacfx[j]=0.0; yvacfy[j]=0.0; yvacfz[j]=0.0;}
  for(k=1;k<=Nz;k++) {zvacf[k]=0.0; zvacfx[k]=0.0; zvacfy[k]=0.0; zvacfz[k]=0.0;}
  denstot=0.0;
  if((t%10)==0) printf("Timestep: %d\n",t);
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     for(l=0;l<=dim-1;l++) 
	        {
	         probv0[i][j][k][l][0]=probv0[i][j][k][l][1]; probv0[i][j][k][l][1]=0.0;
	        }
	     dens[i][j][k][0]=dens[i][j][k][1];
             if(strcmp(densread,"distinguished")==0)
               {
                dens1[i][j][k][0]=dens1[i][j][k][1];
                dens2[i][j][k][0]=dens2[i][j][k][1];
               }
	     ReG[i][j][k][0]=ReG[i][j][k][1];
	     ReG[i][j][k][1]=0.0;
	     ImG[i][j][k][0]=ImG[i][j][k][1];
	     ImG[i][j][k][1]=0.0;
	    }
	}
     }
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     vxmoy=0.0; vymoy=0.0; vzmoy=0.0;
	     if(Lattice[i][j][k]==0)
	       {
	    	/* Neighbouring sites and application of PBC */
	    	i1=i+1;
	    	if(i1>Nx) i1-=Nx;
	    	i2=i-1;
	    	if(i2<=0) i2+=Nx;
	    	j1=j+1;
	    	if(j1>Ny) j1-=Ny;
	    	j2=j-1;
	    	if(j2<=0) j2+=Ny;
	    	k1=k+1;
	    	if(k1>Nz) k1-=Nz;
	    	k2=k-1;
	    	if(k2<=0) k2+=Nz;

                pVtot[i][j][k]=vol_occ[i][j][k];
                pinpore_c[i][j][k]=inpore_c[i][j][k];

	    	/* Probability to arrive in site i,j at t weighted by the velocity */
	    	/* Probability of arriving from another site */
	    	locpfunc=0.0;
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i1][j][k][1]/(kB*T));
		   /* What arrives */
		   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*ptr; 
		   ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*ptr;
		   ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*ptr;
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i1][j][k][1]/(kB*T));
		   /* What arrives */
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*ptr; 
	           locpfunc+=dens[i2][j][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*ptr; 
	           locpfunc+=dens[i][j1][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	   	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*ptr; 
	           locpfunc+=dens[i][j2][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*ptr; 
	           locpfunc+=dens[i][j][k1][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*ptr; 
	           locpfunc+=dens[i][j][k2][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
		/* Adding to the probability of staying */ 
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	locpfunc+=dens[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	/* Integral until t-1 is now calculated */
	    	ReG[i][j][k][1]/=locpfunc; ImG[i][j][k][1]/=locpfunc;
	    	/* Contribution of site i,j at t (product of complex numbers) */
	    	newRe=ReG[i][j][k][1]*cos(wij[i][j][k]*dtime)-ImG[i][j][k][1]*sin(wij[i][j][k]*dtime);
	    	newIm=ReG[i][j][k][1]*sin(wij[i][j][k]*dtime)+ImG[i][j][k][1]*cos(wij[i][j][k]*dtime);
	    	ReG[i][j][k][1]=newRe;
	    	ImG[i][j][k][1]=newIm;
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t*/
	    	/* The path followed does not matter */
	    	parr=0.0;
                parr1=0.0;
                parr2=0.0;
		ptr=exp(-alpha[i1][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i1][j][k][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i1][j][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i1][j][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i1][j][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i1][j][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i1][j][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i1][j][k])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i1][j][k][0]*exp(-(E1site[i][j][k]-E1site[i1][j][k])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i1][j][k])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i1][j][k][0]*exp(-(E2site[i][j][k]-E2site[i1][j][k])/(kB*T))*ptr;
                     }
                  }
		ptr=exp(-alpha[i2][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k])
                  { 
                   parr+=(1/(dir*1.0))*dens[i2][j][k][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i2][j][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i2][j][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i2][j][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i2][j][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i2][j][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i2][j][k])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i2][j][k][0]*exp(-(E1site[i][j][k]-E1site[i2][j][k])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i2][j][k])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i2][j][k][0]*exp(-(E2site[i][j][k]-E2site[i2][j][k])/(kB*T))*ptr;
                     }
                  }
		ptr=exp(-alpha[i][j1][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i][j1][k][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j1][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j1][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j1][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j1][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j1][k])
                  { 
                   parr+=(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j1][k])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i][j1][k][0]*exp(-(E1site[i][j][k]-E1site[i][j1][k])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i][j1][k])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i][j1][k][0]*exp(-(E2site[i][j][k]-E2site[i][j1][k])/(kB*T))*ptr;
                     }
                  }
		ptr=exp(-alpha[i][j2][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i][j2][k][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j2][k])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j2][k][0];
                     }
                   if(E2site[i][j][k]<=E2site[i][j2][k])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j2][k][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j2][k]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j2][k])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i][j2][k][0]*exp(-(E1site[i][j][k]-E1site[i][j2][k])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i][j2][k])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i][j2][k][0]*exp(-(E2site[i][j][k]-E2site[i][j2][k])/(kB*T))*ptr;
                     }
                  }
		ptr=exp(-alpha[i][j][k1][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i][j][k1][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j][k1])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k1][0]; 
                     }
                   if(E2site[i][j][k]<=E2site[i][j][k1])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k1][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j][k1])
                  { 
                   parr+=(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j][k1])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i][j][k1][0]*exp(-(E1site[i][j][k]-E1site[i][j][k1])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i][j][k1])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i][j][k1][0]*exp(-(E2site[i][j][k]-E2site[i][j][k1])/(kB*T))*ptr;
                     }
                  }
		ptr=exp(-alpha[i][j][k2][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2])
                  { 
                   parr+=(1/(dir*1.0))*dens[i][j][k2][0]*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]<=E1site[i][j][k2])
                     {
                      parr1+=ptr*(1/(dir*1.0))*dens1[i][j][k2][0]; 
                     }
                   if(E2site[i][j][k]<=E2site[i][j][k2])
                     {
                      parr2+=ptr*(1/(dir*1.0))*dens2[i][j][k2][0];
                     }
                  }
	    	if(Esite[i][j][k]>Esite[i][j][k2]) 
                  {
                   parr+=(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
                  }
                if(strcmp(densread,"distinguished")==0)
                  {
                   if(E1site[i][j][k]>E1site[i][j][k2])
                     {
                      parr1+=(1/(dir*1.0))*dens1[i][j][k2][0]*exp(-(E1site[i][j][k]-E1site[i][j][k2])/(kB*T))*ptr;
                     }
                   if(E2site[i][j][k]>E2site[i][j][k2])
                     {
                      parr2+=(1/(dir*1.0))*dens2[i][j][k2][0]*exp(-(E2site[i][j][k]-E2site[i][j][k2])/(kB*T))*ptr;
                     }
                  }
	    	/* Probability of leaving site i,j at t */
	    	/* Same as before but the density might have changed */		
	    	pleave=0.0;
                pleave1=0.0;
                pleave2=0.0;
	    	if(Lattice[i1][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][0]/(kB*T));
		  if(Esite[i1][j][k]<=Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i1][j][k]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i1][j][k]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i1][j][k]>Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i1][j][k]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i1][j][k]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i1][j][k]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i1][j][k]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	     	  }
	    	if(Lattice[i2][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][1]/(kB*T));
		  if(Esite[i2][j][k]<=Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i2][j][k]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i2][j][k]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i2][j][k]>Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i2][j][k]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i2][j][k]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i2][j][k]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i2][j][k]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	      	  }
	    	if(Lattice[i][j1][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][2]/(kB*T));
		  if(Esite[i][j1][k]<=Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j1][k]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i][j1][k]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i][j1][k]>Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j1][k]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j1][k]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i][j1][k]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j1][k]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	      	  }
	    	if(Lattice[i][j2][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][3]/(kB*T));
		  if(Esite[i][j2][k]<=Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j2][k]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i][j2][k]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i][j2][k]>Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j2][k]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j2][k]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i][j2][k]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j2][k]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	      	  }
	    	if(Lattice[i][j][k1]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][4]/(kB*T));
		  if(Esite[i][j][k1]<=Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j][k1]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i][j][k1]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i][j][k1]>Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j][k1]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k1]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i][j][k1]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k1]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	      	  }
	    	if(Lattice[i][j][k2]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][5]/(kB*T));
		  if(Esite[i][j][k2]<=Esite[i][j][k]) 
                    {
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j][k2]<=E1site[i][j][k])
                       {
                        pleave1+=ptr*(1/(dir*1.0))*dens1[i][j][k][0];
                       }
                     if(E2site[i][j][k2]<=E2site[i][j][k])
                       {
                        pleave2+=ptr*(1/(dir*1.0))*dens2[i][j][k][0];
                       }
                    }
		  if(Esite[i][j][k2]>Esite[i][j][k])
                    { 
                     pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr;
                    }
                  if(strcmp(densread,"distinguished")==0)
                    {
                     if(E1site[i][j][k2]>E1site[i][j][k])
                       {
                        pleave1+=(1/(dir*1.0))*dens1[i][j][k][0]*exp(-(E1site[i][j][k2]-E1site[i][j][k])/(kB*T))*ptr;
                       }
                     if(E2site[i][j][k2]>E2site[i][j][k])
                       {
                        pleave2+=(1/(dir*1.0))*dens2[i][j][k][0]*exp(-(E2site[i][j][k2]-E2site[i][j][k])/(kB*T))*ptr;
                       }
                    }
	      	  }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
                if(strcmp(densread,"distinguished")==0)
                  {
                   dens1[i][j][k][1]=dens1[i][j][k][0]+parr1-pleave1;
                   dens2[i][j][k][1]=dens2[i][j][k][0]+parr2-pleave2;
                   if(strcmp(iontype,"anion")==0)
                     {
                      vol_occ[i][j][k]=Vanion*dens[i][j][k][1]*part_funct+Vcation*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                   if(strcmp(iontype,"cation")==0)
                     {
                      vol_occ[i][j][k]=Vcation*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vsolvent*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens[i][j][k][1]*part_funct+dens1[i][j][k][1]*part_funct1)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                   if(strcmp(iontype,"solvent")==0)
                     {
                      vol_occ[i][j][k]=Vsolvent*dens[i][j][k][1]*part_funct+Vanion*dens1[i][j][k][1]*part_funct1+Vcation*dens2[i][j][k][1]*part_funct2;
                      inpore_c[i][j][k]=(((dens1[i][j][k][1]*part_funct1+dens2[i][j][k][1]*part_funct2)*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                     }
                  }
                if(strcmp(densread,"undistinguished")==0)
                  {
                   vol_occ[i][j][k]=Vspec*dens[i][j][k][1]*part_funct;
                   inpore_c[i][j][k]=((dens[i][j][k][1]*part_funct*1000.0/Nav)/(2.0*surface_pore[i][j][k]))*graphene_surf;
                  }
	    	/* Average velocity in site i,j */
	    	if((Lattice[i1][j][k]==0)&&(Esite[i1][j][k]<=Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][0]/(kB*T)); 
                 vxmoy=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i1][j][k]==0)&&(Esite[i1][j][k]>Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][0]/(kB*T)); 
                 vxmoy=(1/(dir*1.0))*v0*exp(-(Esite[i1][j][k]-Esite[i][j][k]))*ptr;
                }
	    	if((Lattice[i2][j][k]==0)&&(Esite[i2][j][k]<=Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][1]/(kB*T)); 
                 vxmoy-=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i2][j][k]==0)&&(Esite[i2][j][k]>Esite[i][j][k])) 
                {
                 ptr=exp(-alpha[i][j][k][1]/(kB*T));
                 vxmoy-=(1/(dir*1.0))*v0*exp(-(Esite[i2][j][k]-Esite[i][j][k]))*ptr;
                }
	    	if((Lattice[i][j1][k]==0)&&(Esite[i][j1][k]<=Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][2]/(kB*T)); 
                 vymoy=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i][j1][k]==0)&&(Esite[i][j1][k]>Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][2]/(kB*T)); 
                 vymoy=(1/(dir*1.0))*v0*exp(-(Esite[i][j1][k]-Esite[i][j][k]))*ptr;
                }
	    	if((Lattice[i][j2][k]==0)&&(Esite[i][j2][k]<=Esite[i][j][k]))
                {
                 ptr=exp(-alpha[i][j][k][3]/(kB*T)); 
                 vymoy-=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i][j2][k]==0)&&(Esite[i][j2][k]>Esite[i][j][k]))  
                {
                 ptr=exp(-alpha[i][j][k][3]/(kB*T)); 
                 vymoy-=(1/(dir*1.0))*v0*exp(-(Esite[i][j2][k]-Esite[i][j][k]))*ptr;
                }
	    	if((Lattice[i][j][k1]==0)&&(Esite[i][j][k1]<=Esite[i][j][k])) 
                {
                 ptr=exp(-alpha[i][j][k][4]/(kB*T));
                 vzmoy=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i][j][k1]==0)&&(Esite[i][j][k1]>Esite[i][j][k])) 
                {
                 ptr=exp(-alpha[i][j][k][4]/(kB*T));
                 vzmoy=(1/(dir*1.0))*v0*exp(-(Esite[i][j][k1]-Esite[i][j][k]))*ptr;
                }
	    	if((Lattice[i][j][k2]==0)&&(Esite[i][j][k2]<=Esite[i][j][k])) 
                {
                 ptr=exp(-alpha[i][j][k][5]/(kB*T));
                 vzmoy-=(1/(dir*1.0))*v0*ptr;
                }
	    	if((Lattice[i][j][k2]==0)&&(Esite[i][j][k2]>Esite[i][j][k])) 
                {
                 ptr=exp(-alpha[i][j][k][5]/(kB*T)); 
                 vzmoy-=(1/(dir*1.0))*v0*exp(-(Esite[i][j][k2]-Esite[i][j][k]))*ptr;
                }
	       }
	      lattvacf+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; lattvacf+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      lattvacf+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      lattvacfx+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; lattvacfy+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      lattvacfz+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      xvacf[i]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      xvacfx[i]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; xvacfy[i]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      xvacfz[i]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      yvacf[j]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      yvacfx[j]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; yvacfy[j]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      yvacfz[j]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      zvacf[k]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      zvacfx[k]+=dens[i][j][k][1]*probv0[i][j][k][0][1]*vxmoy; zvacfy[k]+=dens[i][j][k][1]*probv0[i][j][k][1][1]*vymoy; 
	      zvacfz[k]+=dens[i][j][k][1]*probv0[i][j][k][2][1]*vzmoy; 
	      denstot+=dens[i][j][k][1];
	      ReGtot[t]+=dens[i][j][k][1]*ReG[i][j][k][1];
	      ImGtot[t]+=dens[i][j][k][1]*ImG[i][j][k][1];
	  }
	}
    }

  if(Enbar_mode==1)
    {
     for(i=1;i<=Nx;i++)
        {
         for(j=1;j<=Ny;j++)
            {
             for(k=1;k<=Nz;k++)
                {
                 if(Lattice[i][j][k]==0)
                   {
                    for(l=0;l<=5;l++)
                       {
                        if(alpha[i][j][k][l]!=pow(10,20))
                          {
                           alpha[i][j][k][l]=alpha[i][j][k][l]/(kB*T_ref*(pinpore_c[i][j][k]-log(fitcst)));
                           alpha[i][j][k][l]=alpha[i][j][k][l]*(kB*T_ref*(inpore_c[i][j][k]-log(fitcst)));
                           if(inpore_c[i][j][k]-log(fitcst)<0.0) alpha[i][j][k][l]=0.0;
                          }
                       }
                   }
                }
            }
        }
    }

  fprintf(outdenstot,"%e	%lf\n",t*dtime,denstot);
  ReGtot[t]/=denstot;
  ImGtot[t]/=denstot;
  magnitude=sqrt(ReGtot[t]*ReGtot[t]+ImGtot[t]*ImGtot[t]);
  phase=atan(ImGtot[t]/ReGtot[t]);
  fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",t*dtime,ReGtot[t],ImGtot[t],magnitude,phase);
  if((t%nsample)==0)
    {
     ReGtotsmp[treal]=ReGtot[t]; 
     ImGtotsmp[treal]=ImGtot[t]; 
     fprintf(outGsmp,"%e	%lf 	%lf	%lf	%lf\n",t*dtime,ReGtotsmp[treal],ImGtotsmp[treal],magnitude,phase);
     treal++;
     printf("%d\n",treal);
   } 
  fprintf(out,"%e %e\n",t*dtime,lattvacf);
  fprintf(outx,"%e %e\n",t*dtime,lattvacfx);
  fprintf(outy,"%e %e\n",t*dtime,lattvacfy);
  fprintf(outz,"%e %e\n",t*dtime,lattvacfz);
  sumvacf+=lattvacf*dtime/denstot;
  sumvacfx+=lattvacfx*dtime/denstot;
  sumvacfy+=lattvacfy*dtime/denstot;
  sumvacfz+=lattvacfz*dtime/denstot;
  for(i=1;i<=Nx;i++)
     {xsumvacf[i]+=xvacf[i]*dtime; xsumvacfx[i]+=xvacfx[i]*dtime; xsumvacfy[i]+=xvacfy[i]*dtime; xsumvacfz[i]+=xvacfz[i]*dtime;}
  for(j=1;j<=Ny;j++)
     {ysumvacf[j]+=yvacf[j]*dtime; ysumvacfx[j]+=yvacfx[j]*dtime; ysumvacfy[j]+=yvacfy[j]*dtime; ysumvacfz[j]+=yvacfz[j]*dtime;}
  for(k=1;k<=Nz;k++)
     {zsumvacf[k]+=zvacf[k]*dtime; zsumvacfx[k]+=zvacfx[k]*dtime; zsumvacfy[k]+=zvacfy[k]*dtime; zsumvacfz[k]+=zvacfz[k]*dtime;}

  msdvacf[t]=lattvacf*dtime/denstot;
  msdvacfx[t]=lattvacfx*dtime/denstot;
  msdvacfy[t]=lattvacfy*dtime/denstot;
  msdvacfz[t]=lattvacfz*dtime/denstot;
  fprintf(outDt,"%e      %lf\n",t*dtime,(sumvacf/(dim*1.0))/(a*a/(2.0*dim*dtime)));
 }

 /* Calculation of the autocorrelation functions for the signal and sampled signal */
 /* autocorrel(); */
 /* Obtention of the NMR spectrum */
 fourier_transform();
 
 fclose(outdenstot);
 fclose(outG);
 fclose(outGsmp);

 fclose(out);
 fclose(outx);
 fclose(outy);
 fclose(outz);
 fclose(outDt);
 
 summsd1=0.0; summsd1x=0.0; summsd1y=0.0; summsd1z=0.0;
 summsd2=0.0; summsd2x=0.0; summsd2y=0.0; summsd2z=0.0;
 for(t=0;t<=Nmin;t++)
    {
     if(t==0) {msd=0.0; msdx=0.0; msdy=0.0; msdz=0.0;}
     if(t==1)
       {
        msd=1.0*dtime*msdvacf[0];
        msdx=1.0*dtime*msdvacfx[0];
        msdy=1.0*dtime*msdvacfy[0];
        msdz=1.0*dtime*msdvacfz[0];
       }
     if(t>=2)
       {
        summsd1+=msdvacf[t-1];
        summsd1x+=msdvacfx[t-1];
        summsd1y+=msdvacfy[t-1];
        summsd1z+=msdvacfz[t-1];
        summsd2+=(t-1)*dtime*msdvacf[t-1];
        summsd2x+=(t-1)*dtime*msdvacfx[t-1];
        summsd2y+=(t-1)*dtime*msdvacfy[t-1];
        summsd2z+=(t-1)*dtime*msdvacfz[t-1];
        msd=t*dtime*(msdvacf[0]+2.0*summsd1)-2.0*summsd2;
        msdx=t*dtime*(msdvacfx[0]+2.0*summsd1x)-2.0*summsd2x;
        msdy=t*dtime*(msdvacfy[0]+2.0*summsd1y)-2.0*summsd2y;
        msdz=t*dtime*(msdvacfz[0]+2.0*summsd1z)-2.0*summsd2z;
       }
     if(t==Nmin-1) {msdi=msd; msdxi=msdx; msdyi=msdy; msdzi=msdz;}
     if(t==Nmin) {msdf=msd; msdxf=msdx; msdyf=msdy; msdzf=msdz;}
     fprintf(outmsd,"%lf	%lf	%lf	%lf	%lf\n",t*dtime,msd,msdx,msdy,msdz); 
    }
 Df=(msdf-msdi)/dtime; Df/=2.0*dim;
 Dfx=(msdxf-msdxi)/dtime; Dfx/=2.0;
 Dfy=(msdyf-msdyi)/dtime; Dfy/=2.0;
 Dfz=(msdzf-msdzi)/dtime; Dfz/=2.0;
 printf("MSD Diffusion coefficient = %e\n",Df);
 printf("MSD Diffusion coefficient in x direction = %e\n",Dfx);
 printf("MSD Diffusion coefficient in y direction = %e\n",Dfy);
 printf("MSD Diffusion coefficient in z direction = %e\n",Dfz);
 fclose(outmsd);


 /* Printing diffusion coefficient in the x direction */
 out=fopen("Diffusion_coeff_x.dat","w");
 fprintf(out,"#x - D - Dx - Dy - Dz\n");
 for(i=1;i<=Nx;i++)
    {fprintf(out,"%e	%e	%e	%e	%e\n",i*a,xsumvacf[i]/(dim*1.0),xsumvacfx[i],xsumvacfy[i],xsumvacfz[i]);}
 fclose(out);

 /* Printing diffusion coefficient in the y direction */
 out=fopen("Diffusion_coeff_y.dat","w");
 fprintf(out,"#y - D - Dx - Dy - Dz\n");
 for(j=1;j<=Ny;j++)
    {fprintf(out,"%e	%e	%e	%e	%e\n",j*a,ysumvacf[j]/(dim*1.0),ysumvacfx[j],ysumvacfy[j],ysumvacfz[j]);}
 fclose(out);

 /* Printing diffusion coefficient in the z direction */
 out=fopen("Diffusion_coeff_z.dat","w");
 fprintf(out,"#z - D - Dx - Dy - Dz\n");
 for(k=1;k<=Nz;k++)
    {fprintf(out,"%e	%e	%e	%e	%e\n",k*a,zsumvacf[k]/(dim*1.0),zsumvacfx[k],zsumvacfy[k],zsumvacfz[k]);}
 fclose(out);

 printf("VACF Diffusion coefficient: %e\n",sumvacf/(dim*1.0));
 printf("VACF Diffusion coefficient in direction x: %e\n",sumvacfx);
 printf("VACF Diffusion coefficient in direction y: %e\n",sumvacfy);
 printf("VACF Diffusion coefficient in direction z: %e\n",sumvacfz);
 /* Calculates a diffusion anisotropy */
 Dlargest=sumvacfx;
 if(Dlargest<sumvacfy) Dlargest=sumvacfy;
 if(Dlargest<sumvacfz) Dlargest=sumvacfz;
 Dsmallest=sumvacfx;
 if(Dsmallest>sumvacfy) Dsmallest=sumvacfy;
 if(Dsmallest>sumvacfz) Dsmallest=sumvacfz;
 printf("The diffusion anisotropy is (Dlargest-Dsmallest)/Dlargest: %lf\n",(Dlargest-Dsmallest)/Dlargest);
 /* Calculates also the tortuosity */
 Dzero=a*a/(2.0*dim*dtime);
 printf("Bulk diffusion coefficient: %e\n",Dzero);
 printf("Tortuosity: %e\n",Dzero/(sumvacf/(dim*1.0)));

}


/* Spin echo function, very similar to propagation function */
void spin_echo()
{
 FILE *outG;
 int i,j,k,l,t;
 int i1,i2,j1,j2,k1,k2;
 double parr,pleave,denstot;
 double ptr,locpfunc;
 double magnitude,phase; 
 double newRe,newIm;

 /* Total signal not measured and measured */
 sprintf(nameecho,"Signal_spin_echo-%d.dat",tau);
 outG=fopen(nameecho,"w");
 fprintf(outG,"# Time - Real part - Imaginary part - Magnitude - Phase\n");

 /* locpfunc is the "local partition function */
 
 printf("Tau = %d\n",tau);

 /* For each new tau, put ReG/ImG back at 0 */
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++) 
	    {
	     ReG[i][j][k][0]=0.0;
	     ReG[i][j][k][1]=0.0;
	     ImG[i][j][k][1]=0.0;
	     ImG[i][j][k][1]=0.0;
	    }
	}
    } 
 /* Put ReGtot and ImGtot to 0 */
 for(t=1;t<=Nmin;t++)
    {
     ReGtot[t]=0.0;
     ImGtot[t]=0.0;
    }

 /* At t=0, G(0)=1 */
 ReGtot[0]=1.0;
 ImGtot[0]=0.0;
 magnitude=sqrt(ReGtot[0]*ReGtot[0]+ImGtot[0]*ImGtot[0]);
 phase=atan(ImGtot[0]/ReGtot[0]);
 fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",0.0,ReGtot[0],ImGtot[0],magnitude,phase);

 /* Initiation of the propagated function which is
 the probability of arriving at site (i,j) at time 1 
 weighted by the initial velocity v0 */
 denstot=0.0;
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
		/* Neighbouring sites and application of PBC */
		i1=i+1;
		if(i1>Nx) i1-=Nx;
		i2=i-1;
		if(i2<=0) i2+=Nx;
		j1=j+1;
		if(j1>Ny) j1-=Ny;
		j2=j-1;
		if(j2<=0) j2+=Ny;
		k1=k+1;
		if(k1>Nz) k1-=Nz;
		k2=k-1;
		if(k2<=0) k2+=Nz;
	    	/* Probability of arriving in site i,j at t=1 weighted by the velocity */
	    	locpfunc=0.0; 
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
                   /* what arrives */
	           probv0[i][j][k][0][1]-=dens[i1][j][k][0]*v0*ptr; 
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* what stays */
		   /* if it stays, velocity equal zero */
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	          }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
                   /* what arrives */
	           probv0[i][j][k][0][1]-=dens[i1][j][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* what stays */
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	       	   probv0[i][j][k][0][1]+=dens[i2][j][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i2][j][k][0]*ptr;
	       	   locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	       	   probv0[i][j][k][0][1]+=dens[i2][j][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	       	   probv0[i][j][k][1][1]-=dens[i][j1][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j1][k][0]*ptr;
	       	   locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	       	   probv0[i][j][k][1][1]-=dens[i][j1][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	       	   probv0[i][j][k][1][1]+=dens[i][j2][k][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j2][k][0]*ptr;
	       	   locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	      	 }
	    	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	       	   probv0[i][j][k][1][1]+=dens[i][j2][k][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  }
	    	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	       	   probv0[i][j][k][2][1]-=dens[i][j][k1][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j][k1][0]*ptr;
	       	   locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	      	  }
	    	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	       	   probv0[i][j][k][2][1]-=dens[i][j][k1][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  } 
	    	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	       	   probv0[i][j][k][2][1]+=dens[i][j][k2][0]*v0*ptr; 
	       	   locpfunc+=dens[i][j][k2][0]*ptr;
	       	   locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	      	 }
	    	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	      	  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	       	   probv0[i][j][k][2][1]+=dens[i][j][k2][0]*v0*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	       	   locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
		   locpfunc+=dens[i][j][k][0]*(1-ptr);
	      	  }
	     	locpfunc+=(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2])*dens[i][j][k][0];
	        for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t=1 */
	    	parr=0.0;
		ptr=exp(-alpha[i][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k]) parr+=ptr*(1/(dir*1.0))*dens[i1][j][k][0];
	    	if(Esite[i][j][k]>Esite[i1][j][k]) parr+=ptr*(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T));
		ptr=exp(-alpha[i][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k]) parr+=ptr*(1/(dir*1.0))*dens[i2][j][k][0];
	    	if(Esite[i][j][k]>Esite[i2][j][k]) parr+=ptr*(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T));
		ptr=exp(-alpha[i][j][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) parr+=ptr*(1/(dir*1.0))*dens[i][j1][k][0];
	    	if(Esite[i][j][k]>Esite[i][j1][k]) parr+=ptr*(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T));
		ptr=exp(-alpha[i][j][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k]) parr+=ptr*(1/(dir*1.0))*dens[i][j2][k][0];
	    	if(Esite[i][j][k]>Esite[i][j2][k]) parr+=ptr*(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T));
		ptr=exp(-alpha[i][j][k][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) parr+=ptr*(1/(dir*1.0))*dens[i][j][k1][0];
	    	if(Esite[i][j][k]>Esite[i][j][k1]) parr+=ptr*(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T));
		ptr=exp(-alpha[i][j][k][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2]) parr+=ptr*(1/(dir*1.0))*dens[i][j][k2][0];
	    	if(Esite[i][j][k]>Esite[i][j][k2]) parr+=ptr*(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T));
	    	/* Probability of leaving site i,j at t=1 */
	    	pleave=0.0;
	    	if(Lattice[i1][j][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   if(Esite[i1][j][k]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i1][j][k]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T));
	      	  }
	    	if(Lattice[i2][j][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][1]/(kB*T));
		   if(Esite[i2][j][k]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i2][j][k]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T));
		  }
	    	if(Lattice[i][j1][k]==0)
	          {
		   ptr=exp(-alpha[i][j][k][2]/(kB*T));
		   if(Esite[i][j1][k]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i][j1][k]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T));
	          }
	    	if(Lattice[i][j2][k]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][3]/(kB*T));
		   if(Esite[i][j2][k]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i][j2][k]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T));
	          }
	    	if(Lattice[i][j][k1]==0)
	          {
		   ptr=exp(-alpha[i][j][k][4]/(kB*T));
		   if(Esite[i][j][k1]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i][j][k1]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T));
	          }
	    	if(Lattice[i][j][k2]==0)
	      	  {
		   ptr=exp(-alpha[i][j][k][5]/(kB*T));
		   if(Esite[i][j][k2]<=Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0];
		   if(Esite[i][j][k2]>Esite[i][j][k]) pleave+=ptr*(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T));
	          }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
	    	/* Initiation of the function (complex function) allowing for the calculation of the spectra */
	    	ReG[i][j][k][1]=cos(wij[i][j][k]*dtime);		/* Real part */
	    	ImG[i][j][k][1]=sin(wij[i][j][k]*dtime);		/* Imaginary part */
	    	ReGtot[1]+=dens[i][j][k][1]*ReG[i][j][k][1];
	    	ImGtot[1]+=dens[i][j][k][1]*ImG[i][j][k][1];
	       }
	     denstot+=dens[i][j][k][1];
	    }
	}
    }
 ReGtot[1]/=denstot;
 ImGtot[1]/=denstot;
 magnitude=sqrt(ReGtot[1]*ReGtot[1]+ImGtot[1]*ImGtot[1]);
 phase=atan(ImGtot[1]/ReGtot[1]);
 fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",dtime,ReGtot[1],ImGtot[1],magnitude,phase);

 printf("\n");

 /* Propagation 90-tau*/
 for(t=2;t<=tau;t++)
 {
  denstot=0.0;
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     for(l=0;l<=dim-1;l++) 
	        {
	         probv0[i][j][k][l][0]=probv0[i][j][k][l][1]; probv0[i][j][k][l][1]=0.0;
	        }
	     dens[i][j][k][0]=dens[i][j][k][1];
	     ReG[i][j][k][0]=ReG[i][j][k][1];
	     ReG[i][j][k][1]=0.0;
	     ImG[i][j][k][0]=ImG[i][j][k][1];
	     ImG[i][j][k][1]=0.0;
	    }
	}
     }
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
	    	/* Neighbouring sites and application of PBC */
	    	i1=i+1;
	    	if(i1>Nx) i1-=Nx;
	    	i2=i-1;
	    	if(i2<=0) i2+=Nx;
	    	j1=j+1;
	    	if(j1>Ny) j1-=Ny;
	    	j2=j-1;
	    	if(j2<=0) j2+=Ny;
	    	k1=k+1;
	    	if(k1>Nz) k1-=Nz;
	    	k2=k-1;
	    	if(k2<=0) k2+=Nz;
	    	/* Probability to arrive in site i,j at t weighted by the velocity */
	    	/* Probability of arriving from another site */
	    	locpfunc=0.0;
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
		   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*ptr; 
		   ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*ptr;
		   ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*ptr;
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*ptr; 
	           locpfunc+=dens[i2][j][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*ptr; 
	           locpfunc+=dens[i][j1][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	   	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*ptr; 
	           locpfunc+=dens[i][j2][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*ptr; 
	           locpfunc+=dens[i][j][k1][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*ptr; 
	           locpfunc+=dens[i][j][k2][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
		/* Adding to the probability of staying */ 
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	locpfunc+=dens[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	/* Integral until t-1 is now calculated */
	    	ReG[i][j][k][1]/=locpfunc; ImG[i][j][k][1]/=locpfunc;
	    	/* Contribution of site i,j at t (product of complex numbers) */
	    	newRe=ReG[i][j][k][1]*cos(wij[i][j][k]*dtime)-ImG[i][j][k][1]*sin(wij[i][j][k]*dtime);
	    	newIm=ReG[i][j][k][1]*sin(wij[i][j][k]*dtime)+ImG[i][j][k][1]*cos(wij[i][j][k]*dtime);
	    	ReG[i][j][k][1]=newRe;
	    	ImG[i][j][k][1]=newIm;
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t*/
	    	/* The path followed does not matter */
	    	parr=0.0;
		ptr=exp(-alpha[i][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	    	/* Probability of leaving site i,j at t */
	    	/* Same as before but the density might have changed */		
	    	pleave=0.0; 
	    	if(Lattice[i1][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][0]/(kB*T));
		  if(Esite[i1][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i1][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	     	 }
	    	if(Lattice[i2][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][1]/(kB*T));
		  if(Esite[i2][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i2][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j1][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][2]/(kB*T));
		  if(Esite[i][j1][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j1][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j2][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][3]/(kB*T));
		  if(Esite[i][j2][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j2][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k1]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][4]/(kB*T));
		  if(Esite[i][j][k1]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k1]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k2]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][5]/(kB*T));
		  if(Esite[i][j][k2]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k2]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
	       }
	      denstot+=dens[i][j][k][1];
	      ReGtot[t]+=dens[i][j][k][1]*ReG[i][j][k][1];
	      ImGtot[t]+=dens[i][j][k][1]*ImG[i][j][k][1];
	  }
	}
    }
  ReGtot[t]/=denstot;
  ImGtot[t]/=denstot;
  magnitude=sqrt(ReGtot[t]*ReGtot[t]+ImGtot[t]*ImGtot[t]);
  phase=atan(ImGtot[t]/ReGtot[t]);
  fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",t*dtime,ReGtot[t],ImGtot[t],magnitude,phase);
 }

 /* Apply second pulse - Change the phase: for all sites magnitude*e(i*theta)=magnitude*e(i*-theta) 
    theta is the phase */
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
        {
         for(k=1;k<=Nz;k++)
            {
             magnitude=sqrt(ReG[i][j][k][1]*ReG[i][j][k][1]+ImG[i][j][k][1]*ImG[i][j][k][1]);
             phase=atan(ImG[i][j][k][1]/ReG[i][j][k][1]);
             ReG[i][j][k][1]=magnitude*cos(phase);
             ImG[i][j][k][1]=magnitude*sin(-1.0*phase);
            }
        }
    }


 /* Propagation 180 - tau */
 for(t=tau+1;t<=2*tau;t++)
 {
  denstot=0.0;
  /*if((t%10)==0) printf("Timestep: %d\n",t);*/
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     for(l=0;l<=dim-1;l++) 
	        {
	         probv0[i][j][k][l][0]=probv0[i][j][k][l][1]; probv0[i][j][k][l][1]=0.0;
	        }
	     dens[i][j][k][0]=dens[i][j][k][1];
	     ReG[i][j][k][0]=ReG[i][j][k][1];
	     ReG[i][j][k][1]=0.0;
	     ImG[i][j][k][0]=ImG[i][j][k][1];
	     ImG[i][j][k][1]=0.0;
	    }
	}
     }
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
	    	/* Neighbouring sites and application of PBC */
	    	i1=i+1;
	    	if(i1>Nx) i1-=Nx;
	    	i2=i-1;
	    	if(i2<=0) i2+=Nx;
	    	j1=j+1;
	    	if(j1>Ny) j1-=Ny;
	    	j2=j-1;
	    	if(j2<=0) j2+=Ny;
	    	k1=k+1;
	    	if(k1>Nz) k1-=Nz;
	    	k2=k-1;
	    	if(k2<=0) k2+=Nz;
	    	/* Probability to arrive in site i,j at t weighted by the velocity */
	    	/* Probability of arriving from another site */
	    	locpfunc=0.0;
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
		   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*ptr; 
		   ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*ptr;
		   ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*ptr;
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*ptr; 
	           locpfunc+=dens[i2][j][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*ptr; 
	           locpfunc+=dens[i][j1][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	   	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*ptr; 
	           locpfunc+=dens[i][j2][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*ptr; 
	           locpfunc+=dens[i][j][k1][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*ptr; 
	           locpfunc+=dens[i][j][k2][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
		/* Adding to the probability of staying */ 
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	locpfunc+=dens[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	/* Integral until t-1 is now calculated */
	    	ReG[i][j][k][1]/=locpfunc; ImG[i][j][k][1]/=locpfunc;
	    	/* Contribution of site i,j at t (product of complex numbers) */
	    	newRe=ReG[i][j][k][1]*cos(wij[i][j][k]*dtime)-ImG[i][j][k][1]*sin(wij[i][j][k]*dtime);
	    	newIm=ReG[i][j][k][1]*sin(wij[i][j][k]*dtime)+ImG[i][j][k][1]*cos(wij[i][j][k]*dtime);
	    	ReG[i][j][k][1]=newRe;
	    	ImG[i][j][k][1]=newIm;
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t*/
	    	/* The path followed does not matter */
	    	parr=0.0;
		ptr=exp(-alpha[i][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	    	/* Probability of leaving site i,j at t */
	    	/* Same as before but the density might have changed */		
	    	pleave=0.0; 
	    	if(Lattice[i1][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][0]/(kB*T));
		  if(Esite[i1][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i1][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	     	 }
	    	if(Lattice[i2][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][1]/(kB*T));
		  if(Esite[i2][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i2][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j1][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][2]/(kB*T));
		  if(Esite[i][j1][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j1][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j2][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][3]/(kB*T));
		  if(Esite[i][j2][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j2][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k1]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][4]/(kB*T));
		  if(Esite[i][j][k1]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k1]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k2]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][5]/(kB*T));
		  if(Esite[i][j][k2]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k2]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
	       }
	      denstot+=dens[i][j][k][1];
	      ReGtot[t]+=dens[i][j][k][1]*ReG[i][j][k][1];
	      ImGtot[t]+=dens[i][j][k][1]*ImG[i][j][k][1];
	  }
	}
    }
  ReGtot[t]/=denstot;
  ImGtot[t]/=denstot;
  magnitude=sqrt(ReGtot[t]*ReGtot[t]+ImGtot[t]*ImGtot[t]);
  phase=atan(ImGtot[t]/ReGtot[t]);
  fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",t*dtime,ReGtot[t],ImGtot[t],magnitude,phase);
 }
 
 /* Put ReGtot and ImGtot to 0 */
 for(t=1;t<=Nmin;t++)
    {
     ReGtot[t]=0.0;
     ImGtot[t]=0.0;
     ReGtotsmp[t]=0.0;
     ImGtotsmp[t]=0.0;
    }
 
 /* Propagation 2*tau - end; measure signal */
 for(t=1;t<=Nmin;t++)
 {
  denstot=0.0;
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     for(l=0;l<=dim-1;l++) 
	        {
	         probv0[i][j][k][l][0]=probv0[i][j][k][l][1]; probv0[i][j][k][l][1]=0.0;
	        }
	     dens[i][j][k][0]=dens[i][j][k][1];
	     ReG[i][j][k][0]=ReG[i][j][k][1];
	     ReG[i][j][k][1]=0.0;
	     ImG[i][j][k][0]=ImG[i][j][k][1];
	     ImG[i][j][k][1]=0.0;
	    }
	}
     }
  for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     if(Lattice[i][j][k]==0)
	       {
	    	/* Neighbouring sites and application of PBC */
	    	i1=i+1;
	    	if(i1>Nx) i1-=Nx;
	    	i2=i-1;
	    	if(i2<=0) i2+=Nx;
	    	j1=j+1;
	    	if(j1>Ny) j1-=Ny;
	    	j2=j-1;
	    	if(j2<=0) j2+=Ny;
	    	k1=k+1;
	    	if(k1>Nz) k1-=Nz;
	    	k2=k-1;
	    	if(k2<=0) k2+=Nz;
	    	/* Probability to arrive in site i,j at t weighted by the velocity */
	    	/* Probability of arriving from another site */
	    	locpfunc=0.0;
		/* Equivalent to ptr=exp(-alpha[i1][j][k][1]/(kB*T)); */
	    	if((Esite[i][j][k]<=Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
		   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*ptr; 
		   ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*ptr;
		   ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*ptr;
	           locpfunc+=dens[i1][j][k][0]*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i1][j][k])&&(Lattice[i1][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][0]/(kB*T));
		   /* What arrives */
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i1][j][k][0]*probv0[i1][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i1][j][k][0]*ReG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i1][j][k][0]*ImG[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		   /* What stays */
	    	   for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*ptr; 
	           locpfunc+=dens[i2][j][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	    	if((Esite[i][j][k]>Esite[i2][j][k])&&(Lattice[i2][j][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][1]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i2][j][k][0]*probv0[i2][j][k][l][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i2][j][k][0]*ReG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i2][j][k][0]*ImG[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           locpfunc+=dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	    	if((Esite[i][j][k]<=Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*ptr; 
	           locpfunc+=dens[i][j1][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j1][k])&&(Lattice[i][j1][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][2]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j1][k][0]*probv0[i][j1][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j1][k][0]*ReG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j1][k][0]*ImG[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          } 
	   	if((Esite[i][j][k]<=Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*ptr; 
	           locpfunc+=dens[i][j2][k][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j2][k])&&(Lattice[i][j2][k]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][3]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j2][k][0]*probv0[i][j2][k][l][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j2][k][0]*ReG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j2][k][0]*ImG[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           locpfunc+=dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*ptr; 
	           locpfunc+=dens[i][j][k1][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k1])&&(Lattice[i][j][k1]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][4]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k1][0]*probv0[i][j][k1][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k1][0]*ReG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k1][0]*ImG[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
	   	if((Esite[i][j][k]<=Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
		  {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*ptr; 
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*ptr; 
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*ptr; 
	           locpfunc+=dens[i][j][k2][0]*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
	           locpfunc+=dens[i][j][k][0]*(1-exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr);
		  }
	   	if((Esite[i][j][k]>Esite[i][j][k2])&&(Lattice[i][j][k2]==0)) 
	          {
                   ptr=exp(-alpha[i][j][k][5]/(kB*T));
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k2][0]*probv0[i][j][k2][l][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ReG[i][j][k][1]+=dens[i][j][k2][0]*ReG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           ImG[i][j][k][1]+=dens[i][j][k2][0]*ImG[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           locpfunc+=dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	           for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(1-ptr);
	           ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(1-ptr);
	           ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(1-ptr);
	           locpfunc+=dens[i][j][k][0]*(1-ptr);
	          }
		/* Adding to the probability of staying */ 
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]+=dens[i][j][k][0]*probv0[i][j][k][l][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ReG[i][j][k][1]+=dens[i][j][k][0]*ReG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]);
	    	ImG[i][j][k][1]+=dens[i][j][k][0]*ImG[i][j][k][0]*(Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	locpfunc+=dens[i][j][k][0]*(Lattice[i1][j][k]+Lattice[i2][j][k]+Lattice[i][j1][k]+Lattice[i][j2][k]+Lattice[i][j][k1]+Lattice[i][j][k2]);
	    	/* Integral until t-1 is now calculated */
	    	ReG[i][j][k][1]/=locpfunc; ImG[i][j][k][1]/=locpfunc;
	    	/* Contribution of site i,j at t (product of complex numbers) */
	    	newRe=ReG[i][j][k][1]*cos(wij[i][j][k]*dtime)-ImG[i][j][k][1]*sin(wij[i][j][k]*dtime);
	    	newIm=ReG[i][j][k][1]*sin(wij[i][j][k]*dtime)+ImG[i][j][k][1]*cos(wij[i][j][k]*dtime);
	    	ReG[i][j][k][1]=newRe;
	    	ImG[i][j][k][1]=newIm;
	    	for(l=0;l<=dim-1;l++) probv0[i][j][k][l][1]/=locpfunc;
	    	/* Probability of arriving in site i,j at t*/
	    	/* The path followed does not matter */
	    	parr=0.0;
		ptr=exp(-alpha[i][j][k][0]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i1][j][k]) parr+=(1/(dir*1.0))*dens[i1][j][k][0]*exp(-(Esite[i][j][k]-Esite[i1][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][1]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i2][j][k]) parr+=(1/(dir*1.0))*dens[i2][j][k][0]*exp(-(Esite[i][j][k]-Esite[i2][j][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][2]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j1][k]) parr+=(1/(dir*1.0))*dens[i][j1][k][0]*exp(-(Esite[i][j][k]-Esite[i][j1][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][3]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j2][k]) parr+=(1/(dir*1.0))*dens[i][j2][k][0]*exp(-(Esite[i][j][k]-Esite[i][j2][k])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][4]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k1]) parr+=(1/(dir*1.0))*dens[i][j][k1][0]*exp(-(Esite[i][j][k]-Esite[i][j][k1])/(kB*T))*ptr;
		ptr=exp(-alpha[i][j][k][5]/(kB*T));
	    	if(Esite[i][j][k]<=Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*ptr;
	    	if(Esite[i][j][k]>Esite[i][j][k2]) parr+=(1/(dir*1.0))*dens[i][j][k2][0]*exp(-(Esite[i][j][k]-Esite[i][j][k2])/(kB*T))*ptr;
	    	/* Probability of leaving site i,j at t */
	    	/* Same as before but the density might have changed */		
	    	pleave=0.0; 
	    	if(Lattice[i1][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][0]/(kB*T));
		  if(Esite[i1][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i1][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i1][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	     	 }
	    	if(Lattice[i2][j][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][1]/(kB*T));
		  if(Esite[i2][j][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i2][j][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i2][j][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j1][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][2]/(kB*T));
		  if(Esite[i][j1][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j1][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j1][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j2][k]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][3]/(kB*T));
		  if(Esite[i][j2][k]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j2][k]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j2][k]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k1]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][4]/(kB*T));
		  if(Esite[i][j][k1]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k1]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k1]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	if(Lattice[i][j][k2]==0)
	      	 {
                  ptr=exp(-alpha[i][j][k][5]/(kB*T));
		  if(Esite[i][j][k2]<=Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*ptr;
		  if(Esite[i][j][k2]>Esite[i][j][k]) pleave+=(1/(dir*1.0))*dens[i][j][k][0]*exp(-(Esite[i][j][k2]-Esite[i][j][k])/(kB*T))*ptr;
	      	 }
	    	/* Evolution of the density */
	    	/* Initial density plus what arrives - what leaves */
	    	dens[i][j][k][1]=dens[i][j][k][0]+parr-pleave;
	       }
	      denstot+=dens[i][j][k][1];
	      ReGtot[t]+=dens[i][j][k][1]*ReG[i][j][k][1];
	      ImGtot[t]+=dens[i][j][k][1]*ImG[i][j][k][1];
	      ReGtotsmp[t]+=dens[i][j][k][1]*ReG[i][j][k][1];
	      ImGtotsmp[t]+=dens[i][j][k][1]*ImG[i][j][k][1];
	  }
	}
    }
  ReGtot[t]/=denstot;
  ImGtot[t]/=denstot;
  ReGtotsmp[t]/=denstot;
  ImGtotsmp[t]/=denstot;
  magnitude=sqrt(ReGtot[t]*ReGtot[t]+ImGtot[t]*ImGtot[t]);
  phase=atan(ImGtot[t]/ReGtot[t]);
  fprintf(outG,"%e	%lf 	%lf	%lf	%lf\n",(t+2*tau)*dtime,ReGtot[t],ImGtot[t],magnitude,phase);
 }
 
 /* Obtention of the NMR spectrum */
 fourier_transform();

 fclose(outG);

}

void random_lattice()
{
 int i,j,k,n,Nobs,Nfree,l,lim;

 Nobs=Nl*coverage;
 Nfree=Nl-Nobs;
 printf("Nobs %d\n",Nobs);
 lim=10000000;

 if(coverage<=0.5)
 {
 for(n=1;n<=Nobs;n++)
    {
     i=rand()%Nx+1;
     j=rand()%Ny+1;
     k=rand()%Nz+1;
     while((Lattice[i][j][k]==1)&&(l<=lim))
	  {
           j=rand()%Ny+1;
	   l++;
	  } 
     Lattice[i][j][k]=1;
    }
 }
 if(coverage>0.5)
 {
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
        {
         for(k=1;k<=Nz;k++)
            {
             Lattice[i][j][k]=1;
            }
        }
     }
 for(n=1;n<=Nfree;n++)
    {
     i=rand()%Nx+1;
     j=rand()%Ny+1;
     k=rand()%Nz+1;
     while((Lattice[i][j][k]==0)&&(l<=lim))
          {
           j=rand()%Ny+1;
           l++;
          } 
     Lattice[i][j][k]=0;
    }
 }
}


double generlognorm(double mean, double std)
{
 double r1,r2,gauss,lognorm;
 
 r1=rand()%10000;
 r2=rand()%10000;
 r1/=10000;
 r2/=10000;
 /* Gaussian distribution with mean = 0.0 and std = 1.0 */
 gauss=0.0+1.0*sqrt(-2.0*log(1-r1))*cos(2*PI*r2);
 /* Number following a log normal distribution */ 
 lognorm=exp(mean+std*gauss);

 return lognorm;
}

double genergauss(double mean, double std)
{
 double r1,r2,gauss;
 
 r1=rand()%10000;
 r2=rand()%10000;
 r1/=10000;
 r2/=10000;
 gauss=mean+std*sqrt(-2.0*log(1-r1))*cos(2*PI*r2);

 return gauss;

}


/* This function is computed according to Zhang et al. J. Mat. Res. 269, 196-202, 2016 */
/* As for now, only the short-diffusion time regime is described */
double energybarrier(double porsize, double temp, double lattpar, int d)
{

  double e1,e2,e3;
  double Ea;

  e1=4.0/(9.0*sqrt(PI));
  e2=1.0/(1.0*porsize);
  e3=1.0*sqrt(lattpar*lattpar/(2.0*d));
  Ea=-1.0*kB*temp*log(1.0-(e1*e2*e3));

  return Ea;
}


void fourier_transform()
{
 int k,n,count;
 double ReFT,ImFT,*NMR_spec,*NMR_freq;
 double magnitude,phase,freq,save,minvalue,maxvalue,firstfreq,secondfreq;
 FILE *outFT;
 char name[50];

 NMR_spec=dvector(Nvalues);
 NMR_freq=dvector(Nvalues);

 /* Obtention of the NMR spectrum */

 printf("Doing FT.\n");

 sprintf(name,"FT_signal-%d.dat",nsample);
 outFT=fopen(name,"w");
 fprintf(outFT,"# Frequency - Magnitude - Phase - Real part - Imagniary part\n");

 /* Sorting the values according to frequency at the same time as calculating */
 count=0;
 for(k=(Nvalues/2)+1;k<=Nvalues-1;k++)
    {
     ReFT=0.0;
     ImFT=0.0;
     for(n=0;n<=Nvalues-1;n++)
        {
         ReFT+=ReGtotsmp[n]*cos(-2.0*PI*k*n/(Nvalues*1.0));
         ReFT-=ImGtotsmp[n]*sin(-2.0*PI*k*n/(Nvalues*1.0));
         ImFT+=ImGtotsmp[n]*cos(-2.0*PI*k*n/(Nvalues*1.0));
         ImFT+=ReGtotsmp[n]*sin(-2.0*PI*k*n/(Nvalues*1.0));
        }
     freq=-(Nvalues-k)*(1.0/(Nvalues*dwell));
     magnitude=sqrt(ReFT*ReFT+ImFT*ImFT);
     phase=atan(ImFT/ReFT);
     fprintf(outFT,"%lf %lf     %lf     %lf     %lf\n",freq,magnitude/(Nvalues*1.0),phase,ReFT,ImFT);
     NMR_spec[count]=magnitude;
     NMR_freq[count]=freq;
     count++;
    }
 for(k=0;k<=(Nvalues/2);k++)
    {
     ReFT=0.0;
     ImFT=0.0;
     for(n=0;n<=Nvalues-1;n++)
        {
         ReFT+=ReGtotsmp[n]*cos(-2.0*PI*k*n/(Nvalues*1.0));
         ReFT-=ImGtotsmp[n]*sin(-2.0*PI*k*n/(Nvalues*1.0));
         ImFT+=ImGtotsmp[n]*cos(-2.0*PI*k*n/(Nvalues*1.0));
         ImFT+=ReGtotsmp[n]*sin(-2.0*PI*k*n/(Nvalues*1.0));
        }
     freq=k*(1.0/(Nvalues*dwell));
     magnitude=sqrt(ReFT*ReFT+ImFT*ImFT);
     phase=atan(ImFT/ReFT);
     fprintf(outFT,"%lf %lf     %lf     %lf     %lf\n",freq,magnitude/(Nvalues*1.0),phase,ReFT,ImFT);
     NMR_spec[count]=magnitude;
     NMR_freq[count]=freq;
     count++;
    }
 fclose(outFT);

 /* Normalisation of the spectrum */
 /* Finding the minimum */
 minvalue=NMR_spec[0];
 for(k=1;k<=Nvalues-1;k++)
    {
     if(NMR_spec[k]<=minvalue) minvalue=NMR_spec[k];
    }
 /* Removing the baseline */
 for(k=0;k<=Nvalues-1;k++)
    {
     NMR_spec[k]-=minvalue;
    }
 /* integrating */
 save=0.0;
 for(k=0;k<=Nvalues-2;k++)
    {
     /* integral */
     save+=(NMR_spec[k]+NMR_spec[k+1])*(NMR_freq[k+1]-NMR_freq[k])/2.0;
    }
 sprintf(name,"Normalised_spectrum-%d.dat",nsample);
 outFT=fopen(name,"w");
 fprintf(outFT,"# Frequency (ppm) - Magnitude\n");
 maxvalue=0.0;
 for(k=0;k<=Nvalues-1;k++)
    {
     NMR_spec[k]*=larmorfreq/save;
     /* Find the maximum value */
     if(maxvalue<=NMR_spec[k]) maxvalue=NMR_spec[k];
     /* Put frequency in ppm */
     NMR_freq[k]=NMR_freq[k]/larmorfreq;
     fprintf(outFT,"%lf %lf\n",NMR_freq[k],NMR_spec[k]);
    }
 fclose(outFT);
 /* Find the first crossing with half value */
 k=0;
 while(NMR_spec[k]<maxvalue/2.0) k++;
 firstfreq=(NMR_freq[k]+NMR_freq[k-1])/2.0;
 /* Find the first crossing with half value */
 k=Nvalues;
 while(NMR_spec[k]<maxvalue/2.0) k--;
 secondfreq=(NMR_freq[k]+NMR_freq[k+1])/2.0;

 if(tau==0) 
   {
    int_tau0=save;
    printf("FT done (initial integral = %lf).\n",int_tau0);
   }
 if(tau>0) 
   {
    int_tau=save;
    printf("FT done (initial integral = %lf).\n",int_tau);
   }
 printf("Halfwdith (ppm): %lf.\n",secondfreq-firstfreq);
 printf("Halfwdith (Hz): %lf.\n",(secondfreq-firstfreq)*larmorfreq);
 
}

void autocorrel()
{
 FILE *out,*outsmp;
 int i,j,ii,ncorr;
 double Recf,Imcf,Restore,Imstore,magnitude,phase;
 double *norm;

 norm=dvector(Nmin+1);

 /* Need to make it cumulative !!!! */

 out=fopen("Signal_acf.dat","w");
 outsmp=fopen("Signal_acf_smp.dat","w");
 
 fprintf(out,"# Time - Magnitude - Phase - Re - Im\n");
 fprintf(outsmp,"# Time - Magnitude - Phase - Re - Im\n");

 ncorr=Nmin/2;
 for(i=0;i<=Nmin;i++) norm[i]=0.0;

 for(i=1;i<=Nmin;i++)
    {
     if(i<=ncorr)
       {
	Restore=ReGtot[i];
	Imstore=ImGtot[i];
	for(j=i;j>=1;j--)
	   {
	    Recf=Restore*ReGtot[j]-Imstore*ImGtot[j];
	    Imcf=Restore*ImGtot[j]+Imstore*ReGtot[j];
	    ReGcf[i-j]+=Recf;
	    ImGcf[i-j]+=Imcf;
	    norm[i-j]++;
	   }
       }
     if(i<=ncorr)
       {
	ii=i%ncorr;
	Restore=ReGtot[ii];
	Imstore=ImGtot[ii];
	for(j=ncorr;j>=ii+1;j--)
	   {
	    Recf=Restore*ReGtot[j]-Imstore*ImGtot[j];
	    Imcf=Restore*ImGtot[j]+Imstore*ReGtot[j];
	    ReGcf[ii-j+ncorr]+=Recf;
	    ImGcf[ii-j+ncorr]+=Imcf;
	    norm[ii-j+ncorr]++;
	   }
	for(j=ii;j>=1;j--)
	   {
	    Recf=Restore*ReGtot[j]-Imstore*ImGtot[j];
	    Imcf=Restore*ImGtot[j]+Imstore*ReGtot[j];
	    ReGcf[ii-j]+=Recf;
	    ImGcf[ii-j]+=Imcf;
	    norm[ii-j]++;
	   }
       }
    }
 for(i=0;i<=ncorr;i++)
    {
     ReGcf[i]/=(norm[i]*1.0);
     ImGcf[i]/=(norm[i]*1.0);
     magnitude=sqrt(ReGcf[i]*ReGcf[i]+ImGcf[i]*ImGcf[i]);
     phase=atan(ImGcf[i]/ReGcf[i]);
     fprintf(out,"%e	%lf	%lf	%lf	%lf\n",i*dtime,magnitude,phase,ReGcf[i],ImGcf[i]);
    }

 ncorr=Nvalues/2;
 for(i=0;i<=Nmin;i++) norm[i]=0.0;
     
 fclose(out);
 fclose(outsmp);
}


/* Recursive function to find all possible rings for a given atom */
void find_rings(int length)
{
 int i,save_last,j,k,shortcut;

  /* First iteration and if we go bak to the beginning */ 
  if(length==0)
    {
     /*printf("Length 0: exit.\n");*/
     return;
    }
  if(length==1)
   {
    /*printf("Length 1\n");*/
    i=1;
    while((i<=natoms)&&(Neigh_modif[i][current_path[length]]==0)) {i++;}
    /* If there is no neighbouring possibility */
    if(i==(natoms+1))
      {
       length=0; 
       find_rings(length);
      }
    /* Increase the path */
    if(i<(natoms+1))
      {
       length++;
       current_path[length]=i;
       /*printf("length %d, %d %d\n",length,current_path[1],current_path[length]);*/
       /*printf("%d\n",Neigh_modif[6][10]);*/
       Neigh_modif[current_path[length]][current_path[length-1]]=0;
       Neigh_modif[current_path[length-1]][current_path[length]]=0;
       find_rings(length); 
      }
   } 
 /* A solution is found */
 if((current_path[1]==current_path[length])&&(length>3))
   {
    Nrings++;
    printf("New ring found.\n");
    for(j=1;j<=length;j++) printf("%d	",current_path[j]);
    printf("\n----\n");
    Histo_rings[length-1]++;
    PATHS[Nrings][0]=length-1;
    for(j=1;j<=length-1;j++) PATHS[Nrings][j]=current_path[j];
    /* Remove momentarily the neighbours more than one bound away from the central atom */
    if(length>=5)
      {
       for(j=3;j<=length-2;j++) 
          {
           for(k=1;k<=natoms;k++)
	      {
               Neigh_modif[current_path[j]][k]=0;
               Neigh_modif[k][current_path[j]]=0;
	      }
          }
       length=2;
       find_rings(length);
      }
    /* The three-membered ring is a special case */
    if(length==4)
      {
       for(k=1;k<=natoms;k++)
	  {
           Neigh_modif[current_path[2]][j]=0;
           Neigh_modif[j][current_path[2]]=0;
	  }
       length=1;
       find_rings(length);
      } 
   }
 /* Test a new vector */
 if((current_path[1]!=current_path[length])&&(length>1))
   {
    /* If the ring is too long */
    if(length==Lmax) 
      {
       Neigh_modif[current_path[length-1]][current_path[length-2]]=0;
       /*Neigh_modif[current_path[length-2]][current_path[length-1]]=0;*/
       /*printf("%d\n",Neigh_modif[6][10]);*/
       length-=2;
       find_rings(length);
      } 
    i=1;
    while((i<=natoms)&&(Neigh_modif[i][current_path[length]]==0)) {i++;}
    if((i==current_path[length-1])) i=current_path[length-1]+1;
    while((i<=natoms)&&(Neigh_modif[i][current_path[length]]==0)) {i++;}
    /* If the there is no other neighbouring possibility */
    if(i==(natoms+1))
      {
       /*printf("No other neighbour.\n");*/
       Neigh_modif[current_path[length]][current_path[length-1]]=0;
       /*Neigh_modif[current_path[length-1]][current_path[length]]=0;*/
       /*printf("%d\n",Neigh_modif[6][10]);*/
       length--;
       find_rings(length);
      }
    /* Increase the path */
    if(i<(natoms+1))
      {
       length++;
       save_last=current_path[length];
       current_path[length]=i;
       /*for(j=1;j<=length;j++) printf("%d	",current_path[j]);
       printf("\n");*/
       /* Test if there is a shortcut */
       shortcut=0;
       for(j=2;j<length;j++)
	 {
	  if(current_path[j]==i) shortcut=1;
	 }
      /* There is a shortcut, this ring is not valid */
      if(shortcut==1)
        {
         printf("Shortcut\n");
         /*printf("%d %d\n",current_path[length-1],current_path[length-2]);
         printf("%d %d\n",save_last,current_path[length-1]);
         if(current_path[length+1]<=natoms) printf("%d %d\n",current_path[length+1],save_last);*/
         Neigh_modif[current_path[length-1]][current_path[length-2]]=0;
	 Neigh_modif[save_last][current_path[length-1]]=Neigh[save_last][current_path[length-1]];
	 if(current_path[length+1]<=natoms) 
	   {
	    Neigh_modif[current_path[length+1]][save_last]=Neigh[current_path[length+1]][save_last];
	   }
         /*Neigh_modif[current_path[length-1]][current_path[length]]=0;*/
         /*printf("%d\n",Neigh_modif[6][10]);*/
         length-=2;
         find_rings(length);
        }
      /* If the addition leads to a valid option, continue */
      if(shortcut==0)
        {
         find_rings(length); 
        } 
      } 
   }
} 


/* Finds the neighbouring atoms through a simple distance criterion */
void neighbouring()
{
 int i,j,Nneigh,count,ibin;
 double dist,dx,dy,dz;
 FILE *out,*outtable;

 printf("Finding neighbours.\n");

 count=0;
 for(i=1;i<=natoms;i++)
    {
     for(j=i+1;j<=natoms;j++)
	{
	 dx=X[i]-X[j];
	 dy=Y[i]-Y[j];
	 dz=Z[i]-Z[j];
	 /* Periodic boundary conditions */
	 if(dx>(Lx/2.0)) dx-=Lx;
         if(dx<(-Lx/2.0)) dx+=Lx;
	 if(dy>(Ly/2.0)) dy-=Ly;
         if(dy<(-Ly/2.0)) dy+=Ly;
	 if(dz>(Lz/2.0)) dz-=Lz;
         if(dz<(-Lz/2.0)) dz+=Lz;
	 dist=dx*dx+dy*dy+dz*dz;
	 dist=sqrt(dist);
         count++;
	 ibin=floor(dist*nbins/(dmax*1.0));
	 Histo_dist[ibin]++;
	 if(dist<=Rcut) {Neigh[i][j]=1; Neigh[j][i]=1; Neigh_modif[i][j]=1; Neigh_modif[j][i]=1;}
	}
    }

 out=fopen("Histo_dist.dat","w");
 for(i=1;i<=nbins;i++)
    {
     fprintf(out,"%lf	%lf\n",i*dmax/(nbins*1.0),Histo_dist[i]*100.0/(count*1.0));
    }
 fclose(out);

 out=fopen("Neighbours.dat","w");
 outtable=fopen("Neighbours_table.dat","w");
 for(i=1;i<=natoms;i++)
    {
     Nneigh=0;
     fprintf(outtable,"%d	",i);
     for(j=1;j<=natoms;j++)
	{
	 Nneigh+=Neigh[i][j];
         if(Neigh[i][j]==1) fprintf(outtable,"%d	",j);
	}
     fprintf(out,"%d %d\n",i,Nneigh);
     fprintf(outtable,"\n");
    }
 fclose(out);
 fclose(outtable);

 printf("Neigbhours found.\n");

}


/* Calculates the ring statistics */
void ring_stats()
{
 FILE *out,*outpaths;
 int n,i,j,k,count;
 double dx,dy,dz,x0,y0,z0;

 for(n=1;n<=natoms;n++)
    {
     printf("Atom %d\n",n);
     for(i=0;i<=Lmax;i++) current_path[i]=0;
     current_path[1]=n;
     find_rings(1);
     for(i=1;i<=natoms;i++)
	{
	 for(j=i;j<=natoms;j++)
	    {
	     if((i!=n)&&(j!=n))
	       {
	        Neigh_modif[i][j]=Neigh[i][j];
	        Neigh_modif[j][i]=Neigh[j][i];
	       }
	     if((i==n)||(j==n))
	       {
	        Neigh_modif[i][j]=0;
	        Neigh_modif[j][i]=0;
	        Neigh[i][j]=0;
	        Neigh[j][i]=0;
	       }
	    }
	}
    }

 /* Writing results and finding centers of masses of the rings */
 out=fopen("Rings_com.xyz","w");
 outpaths=fopen("Rings_paths.dat","w");
 fprintf(outpaths,"%d\n",Nrings);
 fprintf(out,"%d\n",Nrings);
 fprintf(out,"Step 1\n");
 xcom=dvector(Nrings+1);
 ycom=dvector(Nrings+1);
 zcom=dvector(Nrings+1);
 for(i=1;i<=Nrings;i++)
    {
     k=PATHS[i][0];
     fprintf(outpaths,"%d	",k);
     for(j=1;j<=k;j++)
	{
         if(j==1)
	   {
	    x0=X[PATHS[i][j]];
	    y0=Y[PATHS[i][j]];
	    z0=Z[PATHS[i][j]];
	    xcom[i]=x0*k;
	    ycom[i]=y0*k;
	    zcom[i]=z0*k;
	   }
         if(j>1)
	   {
	    dx=X[PATHS[i][j]]-x0;
	    dy=Y[PATHS[i][j]]-y0;
	    dz=Z[PATHS[i][j]]-z0;
	    /* Periodic boundary conditions */
	    if(dx>(Lx/2.0)) dx-=Lx;
            if(dx<(-Lx/2.0)) dx+=Lx;
	    if(dy>(Ly/2.0)) dy-=Ly;
            if(dy<(-Ly/2.0)) dy+=Ly;
	    if(dz>(Lz/2.0)) dz-=Lz;
            if(dz<(-Lz/2.0)) dz+=Lz;
	    xcom[i]+=dx;
	    ycom[i]+=dy;
	    zcom[i]+=dz;
	   }
	 printf("%d	",PATHS[i][j]);
         fprintf(outpaths,"%d	",PATHS[i][j]);
	}
     xcom[i]/=(k*1.0); ycom[i]/=(k*1.0); zcom[i]/=(k*1.0);
     /* With the method chosen, the mass center can be outside the box at this point. */
     if(xcom[i]>Lx) xcom[i]-=Lx;
     if(xcom[i]<0) xcom[i]+=Lx;
     if(ycom[i]>Lx) ycom[i]-=Ly;
     if(ycom[i]<0) ycom[i]+=Ly;
     if(zcom[i]>Lx) zcom[i]-=Lz;
     if(zcom[i]<0) zcom[i]+=Lz;
     fprintf(out,"F	%lf	%lf	%lf\n",xcom[i],ycom[i],zcom[i]);
     printf("\n");
     fprintf(outpaths,"\n");
    }
 printf("%d rings found.\n",Nrings);
 printf("Done\n");
 fclose(out);
 fclose(outpaths);

 out=fopen("Rings_stats.dat","w");
 for(i=1;i<=Lmax;i++)
    {
     fprintf(out,"%d	%lf\n",i,(Histo_rings[i]*100.0)/(Nrings*1.0));
    }
 fclose(out);

}


/* Function which finds the "normal" vector to each ring, or a vector close to what
   would be the normal (if there is no plane going through all the members of the ring). */
void find_normals()
{
 int i,i2,n,j;
 double A,B,C,D,norm;
 double **histo,angle;
 FILE *out,*outnorm;

 xnorm=dvector(Nrings+1);
 ynorm=dvector(Nrings+1);
 znorm=dvector(Nrings+1);
 /* To calculate the orientations of the rings with respect to ux(1,0,0), uy and uz */
 histo=dmatrix(nzones+1,4);

 for(i=0;i<=nzones;i++)
    {
     for(j=1;j<=3;j++) histo[i][j]=0.0;
    }

 out=fopen("Dir_norm.xyz","w");
 fprintf(out,"%d\n",Nrings);
 fprintf(out,"Step 1\n");
 for(n=1;n<=Nrings;n++)
    {
     /* Find the average plane for the carbon atoms belonging to the rings */
     /* See http://www.les-mathematiques.net/phorum/read.php?13,728203 for the expression */
     A=0; B=0; C=0;
     for(i=1;i<=PATHS[n][0];i++)
	{
	 i2=i+1;
	 if(i2>PATHS[n][0]) i2=1;
	 A+=(Y[PATHS[n][i]]-Y[PATHS[n][i2]])*(Z[PATHS[n][i]]+Z[PATHS[n][i2]]);
	 B+=(Z[PATHS[n][i]]-Z[PATHS[n][i2]])*(X[PATHS[n][i]]+X[PATHS[n][i2]]);
	 C+=(X[PATHS[n][i]]-X[PATHS[n][i2]])*(Y[PATHS[n][i]]+Y[PATHS[n][i2]]);
	}
     D=-A*xcom[n]-B*ycom[n]-C*zcom[n];
     /* Normal to the plane */
     norm=A*A+B*B+C*C;
     norm=sqrt(norm);
     xnorm[n]=A/norm;
     ynorm[n]=B/norm;
     znorm[n]=C/norm;
     fprintf(out,"N  %lf     %lf     %lf\n",xcom[n]+xnorm[n],ycom[n]+ynorm[n],zcom[n]+znorm[n]);
     /* Histogram for the orientation of rings */
     /* With respect to ux */
     angle=calc_angle(xnorm[n],ynorm[n],znorm[n],1,0,0);
     histo[find_anglezone(angle)][1]++;
     /* With respect to uy */
     angle=calc_angle(xnorm[n],ynorm[n],znorm[n],0,1,0);
     histo[find_anglezone(angle)][2]++;
     /* With respect to uz */
     angle=calc_angle(xnorm[n],ynorm[n],znorm[n],0,0,1);
     histo[find_anglezone(angle)][3]++;
    }
 fclose(out);

 /* Output of the orientation histograms */
 out=fopen("Histo_angle.dat","w");
 outnorm=fopen("Histo_angle_norm.dat","w");
 fprintf(out,"# Angle -  with respect to ux - with respect to uy - with respect to uz \n");
 fprintf(outnorm,"# Angle -  with respect to ux - with respect to uy - with respect to uz \n");
 for(i=0;i<=nzones;i++)
    {
     fprintf(out,"%lf	%lf	%lf	%lf\n",i*90.0/(nzones*1.0),histo[i][1]*100.0/(Nrings*1.0),histo[i][2]*100.0/(Nrings*1.0),histo[i][3]*100.0/(Nrings*1.0));
     if(i!=0)
       {
        fprintf(outnorm,"%lf	%lf	%lf	%lf\n",i*90.0/(nzones*1.0),histo[i][1]*100.0/(Nrings*sin(i*(PI/2)/nzones)),histo[i][2]*100.0/(Nrings*sin(i*(PI/2)/nzones)),histo[i][3]*100.0/(Nrings*sin(i*(PI/2)/nzones)));
       }
     if(i==0) fprintf(outnorm,"%lf    %lf     %lf     %lf\n",i*90.0/(nzones*1.0),0.0,0.0,0.0);
    }
 fclose(out);
 fclose(outnorm);

}


/* Calculate the angle between two vectors */
double calc_angle(double coorda1, double coorda2, double coorda3, double coordb1, double coordb2, double coordb3)
{
  double norm1,norm2,angle,scalar_product,save_cos;

  scalar_product=coorda1*coordb1+coorda2*coordb2+coorda3*coordb3;
  norm1=sqrt(coorda1*coorda1+coorda2*coorda2+coorda3*coorda3);
  norm2=sqrt(coordb1*coordb1+coordb2*coordb2+coordb3*coordb3);
  save_cos=scalar_product/(norm1*norm2);
  save_cos=fabs(save_cos);
  angle=180*acos(scalar_product/(norm1*norm2))/PI;

  if(angle>90) angle=180-angle;

  return angle;
}


/* Find the "angle zone" of an angle to plot histograms */
int find_anglezone(double angle)
{
  int izone;

  izone=0;
  while ((angle>izone*90.0/nzones)&&(izone<nzones)) izone++;

  return izone;
}




/******************************* Dynamic memory allocation *****************************/

int **imatrix(int nl, int nc)
{
  int i;
  int **m;
  m=(int **) malloc(nl*sizeof(int*));
  if (m) { m[0]=(int *) malloc(nl*nc*sizeof(int));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

double *dvector(int n)
{
  return (double *) malloc(n*sizeof(double));
}

double **dmatrix(int nl, int nc)
{
  int i;
  double **m;
  m=(double **) malloc(nl*sizeof(double*));
  if (m) { m[0]=(double *) malloc(nl*nc*sizeof(double));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

int *ivector(int n)
{
  return (int *) malloc(n*sizeof(int));
}

double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE)
{
 double ***m;
 int i,j;

 m = (double ***)malloc(sizeof(double **) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (double **)malloc(sizeof(double *) * Y_SIZE);
  
    for (j = 0; j < Y_SIZE; j++)
       m[i][j] = (double *)malloc(sizeof(double) * Z_SIZE);
 }

 return m;
}

double ****fddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE, int W_SIZE)
{
 double ****m;
 int i,j,k;

 m = (double ****)malloc(sizeof(double ***) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (double ***)malloc(sizeof(double **) * Y_SIZE);

    for (j = 0; j < Y_SIZE; j++)
	{
         m[i][j] = (double **)malloc(sizeof(double *) * Z_SIZE);
         for (k = 0; k < Z_SIZE; k++) m[i][j][k] = (double *)malloc(sizeof(double) * W_SIZE);
	}	
 }

 return m;
}

double *****fiveddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE, int W_SIZE, int V_SIZE)
{
 double *****m;
 int i,j,k,l;

 m = (double *****)malloc(sizeof(double ****) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (double ****)malloc(sizeof(double ***) * Y_SIZE);

    for (j = 0; j < Y_SIZE; j++)
	{
         m[i][j] = (double ***)malloc(sizeof(double **) * Z_SIZE);
         for (k = 0; k < Z_SIZE; k++) 
	     {
	      m[i][j][k] = (double **)malloc(sizeof(double *) * W_SIZE);
	      for(l = 0; l < W_SIZE; l++) m[i][j][k][l] = (double *)malloc(sizeof(double) * V_SIZE);
	     }
	}	
 }

 return m;
}

int ***tdimatrix (int X_SIZE, int Y_SIZE, int Z_SIZE)
{
 int ***m;
 int i,j;

 m = (int ***)malloc(sizeof(int **) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (int **)malloc(sizeof(int *) * Y_SIZE);
  
    for (j = 0; j < Y_SIZE; j++)
       m[i][j] = (int *)malloc(sizeof(int) * Z_SIZE);
 }

 return m;
}

