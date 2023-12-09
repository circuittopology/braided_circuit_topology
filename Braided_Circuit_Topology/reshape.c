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

int main(int argc, char *argv[]){
	FILE* fp_lammps;
	FILE* fp_ppa;
	FILE* fp_ic;
	int bufferLength = 255;
	char buffer[bufferLength];
	const char s[2] = " ";
	int nr_atoms = 0;
	int nr_bonds = 0;
	int nr_files = atoi(argv[1]);
	
	char str2[bufferLength]; 
	//int timeslice = 1500000;
	int timeslice = 20000;
	sprintf(str2, "%d", timeslice);
	char buffer1[120], buffer2[120];
	
	for(int k=1; k<=nr_files; k++){
		fp_ic = fopen("chains_ic.data","r");
		snprintf(buffer1, sizeof(char) * 120, "Trj_%d.lammpstrj",k);
		snprintf(buffer2, sizeof(char) * 120, "ppa_ic_%d.data",k);
		fp_lammps = fopen(buffer1, "r");
		fp_ppa = fopen(buffer2,"w");
		
		
		for(int i=0;i<3;i++){
			fgets(buffer,bufferLength,fp_ic);
			if(i==1){
				char *token = strtok(buffer, s);
				nr_atoms = atoi(token);
				//printf("%d\n",nr_atoms);
			} else if(i==2){
				char *token = strtok(buffer, s);
				nr_bonds = atoi(token);
				//printf("%d\n",nr_bonds);
			}
		}
		
		int nr_polymers = nr_atoms-nr_bonds;
		int nr_beads = nr_atoms/nr_polymers;
		
		if(fp_lammps == NULL ) {
			perror ("Error opening LAMMPS file");
	   		return(-1);
	   	}
	   	if(fp_ppa == NULL ) {
			perror ("Error opening PPA file");
	   		return(-1);
	   	}
		
		while(fgets(buffer, bufferLength, fp_lammps)){
			buffer[strcspn(buffer, "\n")] = 0;
			int result = strcmp(buffer, str2);
			if(result==0){
				//printf("T = %s\n",buffer);
				break;
			}
		}
		fprintf(fp_ppa,"LAMMPS data file via write_data\n\n%d atoms\n%d bonds\n\n1 atom types\n1 bond types\n\n",nr_atoms,nr_atoms-nr_polymers);
		for(int i=0; i<7; i++){//Print next lines 7 as information
			fgets(buffer, bufferLength, fp_lammps);
			//printf("%s\n",buffer);
			if(i==3){// x Box sizes
				buffer[strlen(buffer)-1] = '\0';
	   			fprintf(fp_ppa, "%s xlo xhi\n", buffer);
			} else if(i==4){
				buffer[strlen(buffer)-1] = '\0';
	   			fprintf(fp_ppa, "%s ylo yhi\n", buffer);
			} else if(i==5){
				buffer[strlen(buffer)-1] = '\0';
	   			fprintf(fp_ppa, "%s zlo zhi\n\n", buffer);
			}
		}
		fprintf(fp_ppa,"0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 xy xz yz\n\n");
		fprintf(fp_ppa,"Masses\n\n1 1.0\n\nPair Coeffs\n\n1 1 1\n\nBond Coeffs\n\n1 30 1.5 1 1\n\nAtoms\n\n");
		
		while(fgets(buffer, bufferLength, fp_lammps)){
			fprintf(fp_ppa,"%s",buffer);
		}
		
		fprintf(fp_ppa,"\nBonds\n\n");
		
		int bond_nr = 1;
		for(int i=1; i<=nr_polymers; i++){
			for(int j=1; j<nr_beads; j++){
				fprintf(fp_ppa, "%d %d %d %d\n",bond_nr, 1, bond_nr+(i-1), bond_nr+i);
				bond_nr++;
			}
		}
		
		fclose(fp_lammps);
		fclose(fp_ppa);
		fclose(fp_ic);
		}
		
}












