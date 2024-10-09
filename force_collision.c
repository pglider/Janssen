#include "decla.h"

void force_collision(int vlcontact) {
	double ut, utx, uty;
	double x,y;
	double nx, ny;
	double Fn;
	double Vn, Vs, Vsx, Vsy;
	double vsx, vsy;
	double ftx, fty;
	double delta, rij;

	if(vlcontact>=NMAXcontacts){
		printf("\n [WARNING] TOO MANY CONTACST \n");//Should never happen
		getchar();
	}

	//Compute the difference of position between the grains i and j
	//For periodic boundary conditions
	x = CLPdist(disk[j].x -disk[i].x,WIDTH);
	y = disk[j].y -disk[i].y ;

	//Use these for alternative formulations of the normal force. Such formulations can be useful to make some properties independent on the grain size or mass
	//reff=(disk[i].Ray*disk[j].Ray)/(disk[i].Ray+disk[j].Ray);// Effective radius for the dissipation force [OLD FORMULATION]
	//meff=(disk[i].mass*disk[j].mass)/(disk[i].mass+disk[j].mass);// Effective mass for the dissipation force [OLD FORMULATION]

	rij=sqrt(x*x+y*y);
	delta = disk[i].Ray + disk[j].Ray -rij;

	//Overlap abnormally large
	if(delta>R/2) {
		printf("\n [WARNING]@bcltps : %d : abnormally large overlap\n", bcltps);
		printf("\n  (R) = %.3f part i:%d\tet j:%d\tDx: %f\tDy: %f", delta/R, i, j, x, y);
		printf("\n  i:%d\tx: %f\ty: %f", i, disk[i].x, disk[i].y);
		printf("\n  j:%d\tx: %f\ty: %f", j, disk[j].x, disk[j].y);
		err_count+=1;
		//getchar();
	}

	// Definition of the unit vector normal to the contact
	nx=x/rij;
	ny=y/rij;

	// Compute the normal speed Vn=\vec{v}.\vec{n}, viscoelastic model
	Vn=(disk[j].dx-disk[i].dx)*nx + (disk[j].dy-disk[i].dy)*ny;
	Fn=- ( KN* delta - GAMMA*Vn);

	// Component of the sliding speed (relative speed - normal speed + Varignon formula), then norm
	Vsx= disk[j].dx - disk[i].dx - Vn*nx + ny*(disk[j].dOz*disk[j].Ray + disk[i].dOz*disk[i].Ray);
	Vsy= disk[j].dy - disk[i].dy - Vn*ny - nx*(disk[j].dOz*disk[j].Ray + disk[i].dOz*disk[i].Ray);
	Vs= sqrt(Vsx*Vsx+Vsy*Vsy);

	// Cundall's model u_t=\int{\vec{t}.\vec{dl}}
	utx=disk[i].utijx[vlcontact]+KT*Vsx*DT;
	uty=disk[i].utijy[vlcontact]+KT*Vsy*DT;
	ut=sqrt(utx*utx+uty*uty);
	if(fabs(ut/(MU*Fn))>1. && (utx*Vsx)>0 && (uty*Vsy)>0) {
		//First we determine the direction of the normalized sliding speed
		if(Vs==0.) {
			vsx=0.;
			vsy=0.;
		}
		else {
			vsx=Vsx/Vs;
			vsy=Vsy/Vs;
		}
		// Then we determine the friction force
		// The friction force is proportional to the normal force with a coefficient MU and is opposed to the sliding speed (Fn<0 if there is a collision, see above)
		ftx=- MU*Fn*vsx;
		fty=- MU*Fn*vsy;

		disk[i].utijx[vlcontact]=ftx;
		disk[i].utijy[vlcontact]=fty;
	}
	else{
		disk[i].utijx[vlcontact]+=KT*Vsx*DT;
		disk[i].utijy[vlcontact]+=KT*Vsy*DT;
		ftx=disk[i].utijx[vlcontact];//+ 0.001*Vsx;
		fty=disk[i].utijy[vlcontact];//+ 0.001*Vsy;
	}

	//Debugging observed cohesion : Taking fn negative only
	//if(Fn>0){
    //    Fn=0;
	//}

	//Forces and moments
	disk[i].fx+=Fn*nx+ftx;
	disk[i].fy+=Fn*ny+fty;
	disk[i].Mz+= disk[i].Ray * (nx*fty-ny*ftx);

    //Rolling friction
	//disk[i].utijOz[vlcontact]+=(disk[i].dOz-disk[j].dOz)*DT;

	disk[j].fx+=-Fn*nx-ftx;
	disk[j].fy+=-Fn*ny-fty;
	disk[j].Mz+= disk[j].Ray * (nx*fty-ny*ftx);
	disk[j].fnorm+=sqrt((Fn*nx+ftx)*(Fn*nx+ftx) + (Fn*ny+fty)*(Fn*ny+fty));
	
	//[OLD]
	//disk[j].Mz+=CR*(disk[i].dOz-disk[j].dOz);
	//disk[j].Mz-=-KR*(disk[i].utijOz[vlcontact]);

	//Save the forces in a file
	if((DRAWFORCES)&&(bcltps%PASTPS_EPS==0)){
		fprintf(fforce,"%f\t%f\t%f\t%f\t%f\t%d\n",disk[i].x,disk[i].y,disk[j].x,disk[j].y,fabs(Fn),0);
		fprintf(fforce,"%f\t%f\t%f\t%f\t%f\t%d\n",disk[i].x+nx*disk[i].Ray,disk[i].y+ny*disk[i].Ray,disk[i].x+nx*disk[i].Ray+fmin(5*R,R*0.1*ftx/MASS),disk[i].y+ny*disk[i].Ray+fmin(5*R,R*0.1*fty/MASS),fabs(Fn),1);
	}

}

