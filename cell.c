#include "decla.h"

//Returns the grid containing the numbers of particles belonging in each cell, as well as the list of particles being in the same cell

void cell(void) {
	int iicell, jjcell, IXcell, JYcell;

	//Initialize HoC
	for(iicell=0;iicell<=NCELLX;iicell++) {
		for(jjcell=0;jjcell<=NCELLY;jjcell++) {
			HoC[iicell][jjcell]=0;
		}
	}

	for(i=1;i<=N_PART;i++) {
		//Locate particle in cell
        IXcell=(int)(disk[i].x/(grid*R));
        JYcell=(int)(disk[i].y/(grid*R));

        list_link[i]=HoC[IXcell][JYcell];
        HoC[IXcell][JYcell]=i;

        if(IXcell<0 || IXcell>NCELLX){
            printf("[WARNING] X cell error %d i=%d v=%f x=%f bcltps=%d ! \n",IXcell,i,disk[i].dx,disk[i].x,bcltps);
            err_count+=1;
        }
        if(JYcell<0 || JYcell>NCELLY){
            printf("[WARNING] Y cell error %d i=%d v=%f y=%f bcltps=%d ! \n",JYcell,i,disk[i].dy,disk[i].y,bcltps);
            err_count+=1;
        }
	}
}
