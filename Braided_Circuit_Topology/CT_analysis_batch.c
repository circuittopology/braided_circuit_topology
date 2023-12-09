#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>


int contact(int i, int self_cutoff, double cutoff, int nr_atoms, double data[nr_atoms][6]);
void filter_contacts(int nr_atoms, int contactlist[nr_atoms], int filteredlist[nr_atoms], double distancelist[nr_atoms]);
struct information find_mutual_motifs(int i, int j, int nr_atoms, int contactlist[nr_atoms], double data[nr_atoms][6]);
struct motifs find_motifs(int nr_atoms, int contactlist[nr_atoms], double data[nr_atoms][6]);
void swap(int *xp, int *yp);
void delete_duplicates(int number, int a[number]);
void comb(int m, int n, unsigned char *c,int data_array[][13]);
void CT_content(struct information mf, int data_array[][13], int combinations);
int nchoosek(int n,int k);

struct motifs{
	int S, P, X, I2, T2, L2, XL2, I3, T3, I4;
};

struct information{
	int c1, c2, c3, c4, p1, p2, p3, p4; char M[2];
};

int main(int argc, char *argv[]){
	//int nr_polymers = atoi(argv[1]); // # Polymers
	//int nr_beads = atoi(argv[2]); // # Atoms per polymer
	//int nr_atoms = nr_beads*nr_polymers;
	double cutoff_start = strtod(argv[1],NULL);
	double cutoff_end = strtod(argv[2],NULL);
	int RUNS = atoi(argv[3]); // # configurations to analyse
	double kappa = strtod(argv[4],NULL);
	
	FILE* fp_ic;
	fp_ic = fopen("chains_ic.data","r");
	int bufferLength = 255;
	char bufferIC[bufferLength];
	const char s[2] = " ";
	int nr_atoms = 0;
	int nr_bonds = 0;
	
	for(int i=0;i<3;i++){
			fgets(bufferIC,bufferLength,fp_ic);
			if(i==1){
				char *token = strtok(bufferIC, s);
				nr_atoms = atoi(token);
				//printf("%d\n",nr_atoms);
			} else if(i==2){
				char *token = strtok(bufferIC, s);
				nr_bonds = atoi(token);
				//printf("%d\n",nr_bonds);
			}
		}
		
	int nr_polymers = nr_atoms-nr_bonds;	
	char buffer1[120], buffer0[120];
	
	snprintf(buffer0, sizeof(char) * 120, "CT_content_K%0.1f.txt",kappa);
	FILE *fp_content = fopen(buffer0,"w");
	if(fp_content == NULL ) {
    	perror ("Error opening file");
   		return(-1);
   	}
   	
	int combinations = nchoosek(nr_polymers,4);
	if(argc!=5){
		errx(1, "Expecting 4 arguments!");
	}
	FILE* fp_lammps;
	char buffer[bufferLength];
	
	for(int k=1; k<=RUNS; k++){
		int data_array[1][13];
		memset(data_array,0,1*13*sizeof(int));
		unsigned char buf[100];
		comb(nr_polymers, 4, buf,data_array);
		
		for(double cutoff = cutoff_start; cutoff<=cutoff_end; cutoff+=0.01){
			int self_interaction_cutoff = 3; // No self-interaction within # beads
			
			snprintf(buffer1, sizeof(char) * 120, "Trj_%d.lammpstrj",k);
			fp_lammps = fopen(buffer1, "r");
			if(fp_lammps == NULL ) {
				perror ("Error opening file");
		   		return(-1);
		   	}
			char str2[bufferLength]; 
			int timeslice = 20000;
			sprintf(str2, "%d", timeslice);

			while(fgets(buffer, bufferLength, fp_lammps)){
				buffer[strcspn(buffer, "\n")] = 0;
				int result = strcmp(buffer, str2);
				if(result==0){
					//printf("T = %s\n",buffer);
					break;
				}
			}
			for(int i=0; i<7; i++){//Print next lines 7 as information
				fgets(buffer, bufferLength, fp_lammps);
				//printf("%s\n",buffer);
			}
			
			double data[nr_atoms][6];
			memset(data, 0, nr_atoms*6*sizeof(double));
			
			int count_atoms = 0;
			while(fgets(buffer, bufferLength, fp_lammps)){
				char *token = strtok(buffer, " ");
				char *ptr;
				int i = 0;
				while( token != NULL ){
					double strtoken = strtod(token,&ptr);
			  		data[count_atoms][i] = strtoken;
			  		token = strtok(NULL, " ");
			  		i++;
		   		}
		   		count_atoms++;
			}
			fclose(fp_lammps);
			
			int closest;
			double dx, dy, dz;
			int contactlist[nr_atoms], filteredlist[nr_atoms];
			double distancelist[nr_atoms];
			memset(filteredlist,0,nr_atoms*sizeof(int));
			memset(contactlist,0,nr_atoms*sizeof(int));
			
			for(int i = 0; i<nr_atoms; i++){
				closest =  contact(i, self_interaction_cutoff, cutoff, nr_atoms,data);
				contactlist[i]=closest;
					dx = (data[i][3] - data[closest][3]);
					dy = (data[i][4] - data[closest][4]);
					dz = (data[i][5] - data[closest][5]);
				distancelist[i] = dx*dx + dy*dy + dz*dz;
			}
			filter_contacts(nr_atoms, contactlist, filteredlist, distancelist);
			struct motifs mf = find_motifs(nr_atoms, contactlist, data);
			struct information mf_info;
			mf_info.c1 = 0; mf_info.c2 = 0; mf_info.c3 = 0; mf_info.c4 = 0; mf_info.p1 = 0; mf_info.p2 = 0; mf_info.p3 = 0; mf_info.p4 = 0; strcpy(mf_info.M,"NaN");
			
			for(int i=0; i<(nr_atoms-1); i++){
				if(i!=contactlist[i]){ // Only consider distint OTHER contacts, not the contact with itself
					for(int j=i+1; j<nr_atoms; j++){
						if(j!=contactlist[j]){ // Only consider distint OTHER contacts, not the contact with itself
							mf_info = find_mutual_motifs(i, j, nr_atoms, contactlist, data);
							CT_content(mf_info, data_array, combinations);
						}
					}
				}
			}
		}
		
		for(int i=0; i<combinations; i++){
			fprintf(fp_content,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",data_array[i][4],data_array[i][5],data_array[i][6],data_array[i][7],data_array[i][8],data_array[i][9],data_array[i][10],data_array[i][11],data_array[i][12]);
		}
		fflush(stdout);
		printf("\r %d",k);
	}
	fclose(fp_content);
	printf("\n");
	return 0;
}

int contact(int i, int self_cutoff, double cutoff, int nr_atoms, double data[nr_atoms][6]){
	double distsqr = 1000000;
	double distsqrtmp = distsqr;
	int closest_contact = i;
	for(int j=0; j<nr_atoms; j++){
		if(data[j][1]==data[i][1]){//Beads on same polymer
			if(abs(i-j)<self_cutoff){//Within non-self-interaction zone
				distsqr =1000000;
			} else { // Out of non-self-interaction zone --> Can be considered
				double dx = (data[j][3] - data[i][3]);
				double dy = (data[j][4] - data[i][4]);
				double dz = (data[j][5] - data[i][5]);
				distsqr = abs((dx*dx) + (dy*dy) + (dz*dz));
			}
		} else {// Beads on other polymer: no self-interaction prohibition
			double dx = (data[j][3] - data[i][3]);
			double dy = (data[j][4] - data[i][4]);
			double dz = (data[j][5] - data[i][5]);
			distsqr = abs((dx*dx) + (dy*dy) + (dz*dz));
		}
		if((0<distsqr)&&(distsqr<(cutoff*cutoff))){// Close enough for contact!
			if(distsqr<distsqrtmp){// Closer than other contact --> Replaces it
				closest_contact = j;
				distsqrtmp = distsqr;
			}
		}
	}
	return closest_contact;
}

void filter_contacts(int nr_atoms, int contactlist[nr_atoms], int filteredlist[nr_atoms], double distancelist[nr_atoms]){
	int contactindexlist[nr_atoms];
	int count;
	double distancetmp;
	for(int i=0; i<nr_atoms; i++){// Iterate over all atom indices
		memset(contactindexlist,-1,nr_atoms*sizeof(int));
		count = 0;
		distancetmp = distancelist[i];
		for(int j=0; j<nr_atoms; j++){// Search for occurrences of i
			if(contactlist[j]==i){// Occurrence found
				contactindexlist[count] = j; //Add to index list and increase count
				count++;
			}
		}
	}
}

struct motifs find_motifs(int nr_atoms, int contactlist[nr_atoms], double data[nr_atoms][6]){
	struct motifs mf;
	mf.S=0; mf.P=0; mf.X=0; mf.L2=0; mf.T2=0; mf.XL2=0; mf.I2=0; mf.T3=0; mf.I3=0; mf.I4=0;
	int site1, site2, site3, site4, poly1, poly2, poly3, poly4;
	for(int i=0; i<(nr_atoms-1); i++){
		if(contactlist[i]!=i){// check if first contact is distinct
			for(int j=i+1; j<nr_atoms; j++){// Loop over all other contacts
				if(contactlist[j]!=j){// check if second contact is distinct
					site1 = i; site2 = contactlist[i]; // i<j
					site3 = j; site4 = contactlist[j];
					poly1 = (int) data[site1][1]; poly2 = (int) data[site2][1]; poly3 = (int) data[site3][1]; poly4 = (int) data[site4][1];
					if(poly1==poly2){ //Sites 1 and 2 on same polymer, can then only be single-chain CT motifs, T2, I2, or I3
						if(poly3==poly4){ //Sites 3 and 4 also on same polymer, can be only single-chain CT, or I2
							if(poly3==poly1){//All sites on same polymer: Single-chain CT
								if((site2<site3)&&(site2<site4)&&(site1<site4)){//Series motif
									mf.S++;
								} else if(((site4<site1)&&(site1<site2)&&(site2<site3))||((site4<site2)&&(site2<site1)&&(site1<site3))||((site1<site4)&&(site4<site3)&&(site3<site2))||((site1<site3)&&(site3<site4)&&(site4<site2))){//Parallel motif
									mf.P++;
								} else {//Cross motif
									mf.X++;
								}
								
							} else {//On two different polymers --> I2
								mf.I2++;
							}
						} else {// Sites 3 and 4 on different polymers, can be only T2 or I3
							if((poly3==poly1)||(poly4==poly1)){// One site on shared polymer --> T2
								mf.T2++;
							} else {// Sites 3 and 4 on different polymers --> I3
								mf.I3++;
							}
						}
					} else {//Sites 1 and 2 are on different polymers, can then only be L2, XL2, T2, T3, I3, I4
						if(poly3==poly4){// Sites 3 and 4 on same polymer, can only be T2 or I3
							if((poly3==poly1)||(poly3==poly2)){//Sites 3&4 on shared polymer with one of the other contacts 1/2 --> T2
								mf.T2++;
							} else {//Sites 3 & 4 not on shared polymer with 1/2 --> I3
								mf.I3++;
							}
						} else {// Sites 3 and 4 on different polymers, can only be L2, XL2, T3 or I4
							if((poly3!=poly1)&&(poly3!=poly2)&&(poly4!=poly1)&&(poly4!=poly2)){// All on different chains --> I4
								mf.I4++;
							} else if(((poly1==poly3)&&(poly2!=poly4))||((poly1==poly4)&&(poly2!=poly3))||((poly2==poly3)&&(poly1!=poly4))||((poly2==poly4)&&(poly1!=poly3))){///Only one contact coincides, but not the other -->T3
								mf.T3++;
							} else {//Only motif that remains is loop L2 (LX2 is not addressed here).
								mf.L2++;
							} 
						}
					}
				}
			}
		}
	}
	return mf;
}

struct information find_mutual_motifs(int i, int j, int nr_atoms, int contactlist[nr_atoms], double data[nr_atoms][6]){
	struct information mf;
	mf.c1 = i; mf.c2 = contactlist[i]; mf.c3 = j; mf.c4 = contactlist[j]; mf.p1 = (int) data[mf.c1][1]; mf.p2 = (int) data[mf.c2][1]; mf.p3 = (int) data[mf.c3][1]; mf.p4 = (int) data[mf.c4][1];
	int site1, site2, site3, site4, poly1, poly2, poly3, poly4;
		if(contactlist[i]!=i){// check if first contact is distinct
				if(contactlist[j]!=j){// check if second contact is distinct
					site1 = i; site2 = contactlist[i]; // i<j
					site3 = j; site4 = contactlist[j];
					poly1 = (int) data[site1][1]; poly2 = (int) data[site2][1]; poly3 = (int) data[site3][1]; poly4 = (int) data[site4][1];
					if(poly1==poly2){ //Sites 1 and 2 on same polymer, can then only be single-chain CT motifs, T2, I2, or I3
						if(poly3==poly4){ //Sites 3 and 4 also on same polymer, can be only single-chain CT, or I2
							if(poly3==poly1){//All sites on same polymer: Single-chain CT
								if((site2<site3)&&(site2<site4)&&(site1<site4)){//Series motif
									strcpy(mf.M, "S");
								} else if(((site4<site1)&&(site1<site2)&&(site2<site3))||((site4<site2)&&(site2<site1)&&(site1<site3))||((site1<site4)&&(site4<site3)&&(site3<site2))||((site1<site3)&&(site3<site4)&&(site4<site2))){//Parallel motif
									strcpy(mf.M, "P");
								} else {//Cross motif
									strcpy(mf.M, "X");
								}
								
							} else {//On two different polymers --> I2
								strcpy(mf.M, "I2");
							}
						} else {// Sites 3 and 4 on different polymers, can be only T2 or I3
							if((poly3==poly1)||(poly4==poly1)){// One site on shared polymer --> T2
								strcpy(mf.M, "T2");
							} else {// Sites 3 and 4 on different polymers --> I3
								strcpy(mf.M, "I3");
							}
						}
					} else {//Sites 1 and 2 are on different polymers, can then only be L2, XL2, T2, T3, I3, I4
						if(poly3==poly4){// Sites 3 and 4 on same polymer, can only be T2 or I3
							if((poly3==poly1)||(poly3==poly2)){//Sites 3&4 on shared polymer with one of the other contacts 1/2 --> T2
								strcpy(mf.M, "T2");
							} else {//Sites 3 & 4 not on shared polymer with 1/2 --> I3
								strcpy(mf.M, "I3");
							}
						} else {// Sites 3 and 4 on different polymers, can only be L2, XL2, T3 or I4
							if((poly3!=poly1)&&(poly3!=poly2)&&(poly4!=poly1)&&(poly4!=poly2)){// All on different chains --> I4
								strcpy(mf.M, "I4");
							} else if(((poly1==poly3)&&(poly2!=poly4))||((poly1==poly4)&&(poly2!=poly3))||((poly2==poly3)&&(poly1!=poly4))||((poly2==poly4)&&(poly1!=poly3))){///Only one contact coincides, but not the other -->T3
								strcpy(mf.M, "T3");
							} else {//Only motif that remains is loop L2 (LX2 is not addressed here).
								strcpy(mf.M, "L2");
							} 
						}
					}
				}
		}
	return mf;
}

void comb(int m, int n, unsigned char *c, int data_array[][13]){
	int i;
	int row = 0;
	int tmp;
	for (i = 0; i < n; i++){
		c[i] = n - i;
	}
	while(1){
		for (i = n; i--;){
			tmp = (int) (c[i]);
			data_array[row][n-i-1] = tmp;
		}
		row++;
		/* this check is not strictly necessary, but if m is not close to n,
		   it makes the whole thing quite a bit faster */
		i = 0;
		if (c[i]++ < m) continue;
		for (; c[i] >= m - i;) if (++i >= n) return;
		for (c[i]++; i; i--) c[i-1] = c[i] + 1;
	}
}


void CT_content(struct information mf, int data_array[][13], int combinations){
	char S[] = "S", P[] = "P", X[] = "X";
	char L2[] = "L2", T2[] = "T2", I2[] = "I2";
	char T3[] = "T3", I3[] = "I3", I4[] = "I4";
	
	if((strcmp(mf.M,S)==0)||(strcmp(mf.M,P)==0)||(strcmp(mf.M,X)==0)){// Contact is S/P/X
		for(int i=0; i<combinations; i++){
			if((data_array[i][0]==mf.p1)||(data_array[i][1]==mf.p1)||(data_array[i][2]==mf.p1)||(data_array[i][3]==mf.p1)){ // Check which combinations involve chain p1
				if(strcmp(mf.M,S)==0){ // Add to S
					data_array[i][4]+=1;
				} else if(strcmp(mf.M,P)==0){ // Add to P
					data_array[i][5]+=1;
				} else if(strcmp(mf.M,X)==0){ // Add to X
					data_array[i][6]+=1;
				}
			}
		}
	} else if((strcmp(mf.M,L2)==0)||(strcmp(mf.M,T2)==0)||(strcmp(mf.M,I2)==0)){ //Contact is L2/T2/I2
		int p1 = mf.p1;
		int p2 = mf.p2;
		if(mf.p1!=mf.p2){
			p2 = mf.p2;
		} else if(mf.p1!=mf.p3){
			p2 = mf.p3;
		} else if(mf.p1!=mf.p4){
			p2 = mf.p4;
		}
		if(p1>p2){//Reverse
			int ptmp = p2;
			p2 = p1;
			p1 = ptmp;
		}
		for(int i=0; i<combinations; i++){
			if((data_array[i][0]==p1 && data_array[i][1]==p2)||(data_array[i][0]==p1 && data_array[i][2]==p2)||(data_array[i][0]==p1 && data_array[i][3]==p2)||(data_array[i][1]==p1 && data_array[i][2]==p2)||(data_array[i][1]==p1 && data_array[i][3]==p2)||(data_array[i][2]==p1 && data_array[i][3]==p2)){
				if(strcmp(mf.M,L2)==0){ // Add to L2
					data_array[i][7]+=1;
				} else if(strcmp(mf.M,I2)==0){ // Add to I2
					data_array[i][8]+=1;
				} else if(strcmp(mf.M,T2)==0){ // Add to T2
					data_array[i][9]+=1;
				}
			}
		}
	} else if((strcmp(mf.M,T3)==0)||(strcmp(mf.M,I3)==0)){ //Contact is T3/I3
		int p1 = mf.p1, p2 =mf.p2, p3=mf.p3;
		if(mf.p2==p1){ // 1 and 2 are the same, remaining strands 3 and 4 make up the set --> done
			p2 = mf.p3;
			p3 = mf.p4;
		} else {// 1 and 2 are not the same --> 1 and 2 are part of the set, check other strands
			int p2 = mf.p2;
			if((mf.p3==p2)||(mf.p3==p1)){ // 2 and 3 are the same, strands 1, 2 and 4 make up the set --> done
				p3 = mf.p4;
			} else {
				p3 = mf.p3;
			}
		}
		if (p1 > p3){
   			swap(&p1, &p3);
   		}
		if (p1 > p2){
   			swap(&p1, &p2);
   		}
		if (p2 > p3){
		   swap(&p2, &p3);
		}
		for(int i=0; i<combinations; i++){
			if((data_array[i][0]==p1 && data_array[i][1]==p2 && data_array[i][2]==p3) || (data_array[i][0]==p1 && data_array[i][1]==p2 && data_array[i][3]==p3) || (data_array[i][0]==p1 && data_array[i][2]==p2 && data_array[i][3]==p3) || (data_array[i][1]==p1 && data_array[i][2]==p2 && data_array[i][3]==p3)){
				if(strcmp(mf.M,T3)==0){ // Add to T3
					data_array[i][10]+=1;
				} else if(strcmp(mf.M,I3)==0){ // Add to I3
					data_array[i][11]+=1;
				}
			}
		}
	} else if(strcmp(mf.M,I4)==0){//Contact is I4
		int p1 = mf.p1, p2 = mf.p2, p3 =mf.p3, p4=mf.p4;
		if(p1>p2) swap(&p1,&p2);
		if(p3>p4) swap(&p3,&p4);
		if(p1>p3) swap(&p1,&p3);
		if(p2>p4) swap(&p2,&p4);
		if(p2>p3) swap(&p2,&p3);
		for(int i=0; i<combinations; i++){
			if(data_array[i][0]==p1 && data_array[i][1]==p2 && data_array[i][2]==p3 && data_array[i][3]==p4){
				data_array[i][12]+=1;
			}
		}
	}
}

void swap(int *xp, int *yp){
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

int nchoosek(int n,int k){
	if(k == 0){
		return 1;
	} else {
		return (n*nchoosek(n-1, k-1))/k;
	}
}
