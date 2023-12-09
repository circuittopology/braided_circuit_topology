#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define xl -8.0
#define xh 8.0 // Long direction
#define yl -8.0
#define yh 8.0
#define zl -8.0
#define zh 8.0

#define FENE 0.972

int main(int argc, char *argv[]){
	int n = atoi(argv[1]); // # Polymers
	int m = atoi(argv[2]); // # Atoms per polymer
	int by = atoi(argv[3]); // Initial box size y
	int bz = atoi(argv[4]); // Initial box size z
	
	if(argc!=5){
		errx(1, "Expecting 4 arguments!");
	}
	struct timeval time;
			gettimeofday(&time,NULL);
	const gsl_rng_type *gsl_rng_env_setup();
			gsl_rng *r;
			r=gsl_rng_alloc (gsl_rng_default);
			gsl_rng_set(r,time.tv_usec);
	
	double polymer_length = m*FENE;
	
	if(polymer_length>= abs(xh-xl)){
		errx(1, "Warning: polymer is too long for box. LAMMPS is going to give a FENE error. Adjust box sizes.");
	}
	
	FILE *fp_ic = fopen("chains_ic.data","w");
	if(fp_ic == NULL ){
    	perror ("Error opening file");
   		return(-1);
   	}
	fprintf(fp_ic, "Pure Nylon6\n%d atoms\n%d bonds\n%d angles\n%d dihedrals\n\n", n*m, (m-1)*n, (m-2)*n, 0);
	fprintf(fp_ic, "%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n\n", 1, 1, 1, 0);
	fprintf(fp_ic, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n\nMasses\n\n1 1\n\nAtoms\n\n", xl, xh, yl, yh, zl, zh);
	
	for(int i=1; i<=n; i++){
		double Ry = gsl_ran_flat(r, (double) -by/2.0, (double) by/2.0);
		double Rz = gsl_ran_flat(r, (double) -bz/2.0, (double) bz/2.0);
		double Rxinit = (double) xl;
		for(int j=1; j<=m; j++){
			double Rx = Rxinit + (j-1)*FENE;
			fprintf(fp_ic,"%d %d %d %lf %lf %lf\n", (i-1)*m + j, i, 1, Rx , Ry, Rz);
		}
	}
	fprintf(fp_ic,"\nBonds\n\n");
	
	int bond_nr = 1;
	
	for(int i=1; i<=n; i++){
		for(int j=1; j<m; j++){
			fprintf(fp_ic, "%d %d %d %d\n",bond_nr, 1, bond_nr+(i-1), bond_nr+i);
			bond_nr++;
		}
	}
	fprintf(fp_ic,"\nAngles\n\n");
	
	int angle_nr = 1;
	for(int i=1; i<=n; i++){
		for(int j=1; j<(m-1); j++){
			fprintf(fp_ic, "%d %d %d %d %d\n",angle_nr, 1, angle_nr+(2*(i-1)), angle_nr+(2*(i-1))+1, angle_nr+(2*(i-1))+2);
			angle_nr++;
		}
	}
	fclose(fp_ic);	
}
