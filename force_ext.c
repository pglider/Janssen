#include "decla.h"

void force_ext(void) {
	//Implement gravity for grains awake
	if(disk[i].fixed==0){
        disk[i].fy = disk[i].fy - G*disk[i].mass*cos(theta * PI/180);
        disk[i].fx = disk[i].fx + G*disk[i].mass*sin(theta * PI/180);
	}
	else{
        disk[i].fy=0;
        disk[i].fx=0;
	}

    //Limits of space allowed
    //limbas=0;
    //limd=L*R*grid;
    //limg=0;

    //Bouncing force against walls of allowed space
    //if(disk[i].x + disk[i].Ray > limd && disk[i].y>limbas){
    //        disk[i].fx = disk[i].fx + KN*(limd -disk[i].x - disk[i].Ray);
    //}
    //if(disk[i].x - disk[i].Ray < limg && disk[i].y>limbas)
    //    disk[i].fx = disk[i].fx + KN*(limg -disk[i].x + disk[i].Ray);
}

void force_cohesion(void) {
    double x,y;
    double nx, ny, rij;

    if(i>=p0 && (i-p0)%2==0){
        x = CLPdist(disk[i+1].x -disk[i].x,WIDTH);
        y = disk[i+1].y -disk[i].y ;
        rij=sqrt(x*x+y*y);
        nx=x/rij;
        ny=y/rij;
        disk[i].fx+=nx*FORCE_COHESION;
        disk[i].fy+=ny*FORCE_COHESION;
        disk[i+1].fx-=nx*FORCE_COHESION;
        disk[i+1].fy-=ny*FORCE_COHESION;
    }

}
