#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double kb=1.380e-23;
const int nAtoms=864;
const int nUnit=6;

void populateBox(double (*crd)[3], double (*unitcrd)[3], int nUnit, double unitbox){
	int i = 0;
	for(int nx=0;nx<nUnit;nx++) { //iterate over x coords
		for(int ny=0;ny<nUnit;ny++) { //iterate over y coords
			for(int nz=0;nz<nUnit;nz++) { //iterate over z coords
				for(int j=0;j<4;j++) { //iterate over 4 atoms of FCC
					crd[i][0] = (unitcrd[j][0] + nx)*unitbox;
					crd[i][1] = (unitcrd[j][1] + ny)*unitbox;
					crd[i][2] = (unitcrd[j][2] + nz)*unitbox;
					//printf("%12.4e %12.4e %12.4e\n",crd[i][0], crd[i][1], crd[i][2]);
					i++;
				}
			}
		}
	}	
}

int moveParticle(int nAtoms, double maxMCsize, double* move) {
	double MCsize;
	int sel;
	double length;
	length=2.0;
	MCsize=((double)rand())/RAND_MAX*maxMCsize; //random number between 0 and maxMCsize (3e-10)
	sel=floor(((double)rand())/(RAND_MAX-1)*nAtoms); //random number between 0 and nAtoms (486)
	while(length>1.0) {
		move[0]=2*((double)rand())/RAND_MAX-1;
		move[1]=2*((double)rand())/RAND_MAX-1;
		move[2]=2*((double)rand())/RAND_MAX-1;
		length=move[0]*move[0]+move[1]*move[1]+move[2]*move[2];
	}
	length=sqrt(length);
	move[0]=move[0]/length*MCsize;
	move[1]=move[1]/length*MCsize;
	move[2]=move[2]/length*MCsize;
	return sel;
}

double getDistSq(double* a, double* b, double box) {
	double distSq = 0.0;
	double link[3];
	for(int m=0;m<3;m++) {
		link[m]=a[m]-b[m];
		//pbc
		if(link[m]>box/2) {
			link[m]-=box;
		} else if(link[m]<=-1*box/2) {
			link[m]+=box;
		}
		distSq+=link[m]*link[m];
	}
	return distSq;
}
double getE(int sel, int nAtoms, double* crd, double (*crd2)[3], double box, double epsilon, double sigma6, double sigma12) {
	double distSq, dist6, dist12;
	double link[3];
	double e = 0.0;
	for(int k=0;k<nAtoms;k++) {
		if(k!=sel) {
			distSq = getDistSq(crd, crd2[k], box);
			dist6=distSq*distSq*distSq;
			dist12=dist6*dist6;
			e+=4*epsilon*((sigma12/dist12)-(sigma6/dist6));
		}
	}
	return e;
}

void getTryCrd(double* tryCrd, double* crd, double* move, double box) {
	for(int m=0;m<3;m++) {
		tryCrd[m]=crd[m]+move[m];
		if(tryCrd[m]>box) {
			tryCrd[m]-=box;
		} else if(tryCrd[m]<0) {
			tryCrd[m]+=box;
		}
	}
}

int main() {
	FILE *io;
	double temp=100;
	double crd[nAtoms][3];
	double box=3.46812e-9;
	double unitbox=box/6;	
	double unitcrd[4][3];
	int nx,ny,nz;
	int i,j,k,l,m;
	double sigma=3.4e-10;
	double sigmaSq,sigma6,sigma12;
	double epsilon=kb*120;
	int nMC=10000;
	double maxMCsize=3.0e-10;
	double MCsize;
	double move[3];
	double tryCrd[3];
	double length;
	int sel;
	double eOld,eNew,eDelta;
	double distSq,dist6,dist12;
	double link[3];
	double cut=9.0e-10;
	double prob,dice;
	int trjOut=10;

	sigmaSq=sigma*sigma;
	sigma6=sigmaSq*sigmaSq*sigmaSq;
	sigma12=sigma6*sigma6;

	unitcrd[0][0]=0.0;
	unitcrd[0][1]=0.0;
	unitcrd[0][2]=0.0;
	unitcrd[1][0]=0.5;
    unitcrd[1][1]=0.5;
    unitcrd[1][2]=0.0;
	unitcrd[2][0]=0.5;
    unitcrd[2][1]=0.0;
    unitcrd[2][2]=0.5;
	unitcrd[3][0]=0.0;
    unitcrd[3][1]=0.5;
    unitcrd[3][2]=0.5;
	
	populateBox(crd, unitcrd, nUnit, unitbox);
	io=fopen("MC-w-E.xyz","w");
	for(i=1;i<=nMC;i++) {
		for(j=0;j<nAtoms;j++) {
			sel = moveParticle(nAtoms, maxMCsize, move);
			eOld=0.0;
			eOld = getE(sel, nAtoms, crd[sel], crd, box,  epsilon, sigma6, sigma12);
			getTryCrd(tryCrd, crd[sel], move, box);
			eNew=0.0;
			eNew = getE(sel, nAtoms, tryCrd, crd, box, epsilon, sigma6, sigma12);
			eDelta=eNew-eOld;
			if(eDelta<0.0) {
				for(m=0;m<3;m++) {
					crd[sel][m]=tryCrd[m];
				}
			} else {
				prob=exp(-1*eDelta/(kb*temp));
				dice=((double)rand())/RAND_MAX;
				if(dice<=prob) {
					for(m=0;m<3;m++) {
						crd[sel][m]=tryCrd[m];
					}
				}
			}
		}
		if(i%trjOut==0) {
			fprintf(io,"%d\n\n",nAtoms);
			for(j=0;j<nAtoms;j++) {
				fprintf(io,"X %10.3f %10.3f %10.3f\n",crd[j][0]*1e10,crd[j][1]*1e10,crd[j][2]*1e10);
				//fprintf(io, "%10.3f",eDelta*1e20);
			}
			printf("Cycle %d of %d\n", i, nMC);
			fflush(stdout);
		}
	}
	fclose(io);

	return 0;
}
