#include "decla.h"
// Stocke des informations dans un fichier

void save_data(void){
	// Crée et ouvre un fichier
	int ips = bcltps/POS_FREQ;
	char namepos[255], filepos[255]="";
	FILE *fpos;
	sprintf(namepos,"pos_%.7d.txt", ips);
	strcat(filepos,save_folder);
	strcat(filepos,namepos);
	//printf("filepos=%s\n",filepos);
	fpos=fopen(filepos,"w");

	// L'énergie cinétique est la somme des énergies cinétiques de chacun des grains

	// On stocke les rayons, positions et vitesses des particules dans le fichier posXXXX.txt
	//fprintf(fpos,"%d\t-1\t%lf\t%lf\t%lf\n", bcltps, X_slider, Y_slider, 1000*Fx_slider);
	for(i=1;i<=N_PART;i++) {
		fprintf(fpos,"%d\t%.2e\t%.7lf\t%.7lf\t%.2e\t%.2e\t%.2e\t%d\t%.2e\t%.2e\t%.2e\n", i, disk[i].Ray, disk[i].x, disk[i].y, disk[i].dx, disk[i].dy, disk[i].dOz, disk[i].Nb_Contact, disk[i].fx, disk[i].fy, disk[i].fnorm);
	}
	fclose(fpos);
}
