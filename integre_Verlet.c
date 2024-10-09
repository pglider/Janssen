#include "decla.h"

// Intégration des équations du mouvement par la méthode de Verlet.

void integre_Verlet(void){
	// Integrating the equations of motion for the grains: position and angular velocities
    double xnew, ynew, Oznew;

	//Using Verlet's integration method
	
    //Note : this bit should be useless
    if(disk[i].fixed){
        // xnew=disk[i].x;
        // ynew= disk[i].y;
        // Oznew=disk[i].Oz;
        disk[i].fx=0;
        disk[i].fy=0;
        disk[i].Mz=0;
        disk[i].dx=0;
        disk[i].dy=0;
        disk[i].dOz=0;
        disk[i].xold=disk[i].x;
        disk[i].yold=disk[i].y;
        disk[i].Ozold=disk[i].Oz;
    }
    else {
        xnew=CLPdist(2*disk[i].x-disk[i].xold, WIDTH)+DT*DT*disk[i].fx /disk[i].mass;
        ynew=2*disk[i].y-disk[i].yold+DT*DT*disk[i].fy/disk[i].mass;

        // Oznew=disk[i].Oz;
        Oznew=2*disk[i].Oz-disk[i].Ozold+DT*DT*disk[i].Mz /disk[i].II;

        xnew = CLP(xnew,WIDTH);
        disk[i].x=CLP(disk[i].x,WIDTH);
        disk[i].xold=CLP(disk[i].xold,WIDTH);

    //	disk[i].dx=(xnew-disk[i].xold)/(2*DT);
        disk[i].dx=CLPdist((xnew-disk[i].xold),WIDTH)/(2*DT);
        disk[i].dy=(ynew-disk[i].yold)/(2*DT);

        disk[i].dOz=(Oznew-disk[i].Ozold)/(2*DT);
        // disk[i].dOz=0;

        disk[i].xold=disk[i].x;
        disk[i].yold=disk[i].y;

        disk[i].Ozold=disk[i].Oz;

        disk[i].x=xnew;
        disk[i].y=ynew;

        disk[i].Oz=Oznew;
        //Reset forces for next computation
        disk[i].fx=0;
        disk[i].fy=0;
        disk[i].Mz=0;
    }

}

void integre_VerletStop(void){
    //Do not integrate the equations of motion for the grains with a fixed position

    disk[i].dx=0;
    disk[i].dy=0;

    disk[i].dOz=0;

    disk[i].xold=disk[i].x;
    disk[i].yold=disk[i].y;

    disk[i].Ozold=disk[i].Oz;

    disk[i].Mz=0.;

    disk[i].fx=0.;
    disk[i].fy=0.;
    disk[i].fnorm=0.;
}

