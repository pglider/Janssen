#include "decla.h"

void force_mag(double mag_mu) {
	
	int IXcell, JYcell, ixTemp, jyTemp, ixTempbis;
	double x,y,rij, nx, ny;
	double mix, miy, mjx, mjy;
	double mixn, miyn, mjxn, mjyn;
	double Fx, Fy, torque_j;

	// identify the cell part. i belongs to
	IXcell=(int) (disk[i].x/(grid*R));
	JYcell=(int) (disk[i].y/(grid*R));

	// look for possible interactions with part. in the neighbor cells
	for(ixTemp=IXcell-MAG_NEIGHBOR_CELLS;ixTemp<=IXcell+MAG_NEIGHBOR_CELLS ; ixTemp++){
		for(jyTemp=JYcell-MAG_NEIGHBOR_CELLS;jyTemp<=JYcell+MAG_NEIGHBOR_CELLS ; jyTemp++){
			ixTempbis=CLPcell(ixTemp,NCELLX);

			if(ixTempbis>=0 && jyTemp>=0){

				j=HoC[ixTempbis][jyTemp];// j = Head of Chain
				while(j!=0){// j=0 means end of chain
					if(i!=j){
						// magnetic?
						if(disk[j].ismag==1){  

							Fx = 0;
							Fy = 0;
							torque_j = 0;
							//torque_i = 0;

							// Compute the difference of position between the grains i and j
							//x = disk[j].x -disk[i].x;
							x = CLPdist(disk[j].x -disk[i].x,WIDTH);
							y = disk[j].y -disk[i].y ;

							mixn = cos(disk[i].Oz);
							miyn = sin(disk[i].Oz);
							mjxn = cos(disk[j].Oz);
							mjyn = sin(disk[j].Oz);

							mix = mag_mu * mixn;
							miy = mag_mu * miyn;
							mjx = mag_mu * mjxn;
							mjy = mag_mu * mjyn;

							rij=sqrt(x*x+y*y);// Algebraic distance between the grains
							// Definition of the unit vector normal to the interaction
							nx=x/rij;
							ny=y/rij;

							Fx = PREF1/(rij*rij*rij*rij)*((ny*mix - nx*miy) * mjy + (ny*mjx - nx*mjy) * miy - 2*nx*(mix*mjx + miy*mjy) + 5*nx*(ny*mix - nx*miy)*(ny*mjx - nx*mjy));
							Fy = PREF1/(rij*rij*rij*rij)*((nx*miy - ny*mix) * mjx + (nx*mjy - ny*mjx) * mix - 2*ny*(mix*mjx + miy*mjy) + 5*ny*(ny*mix - nx*miy)*(ny*mjx - nx*mjy));

							//torque_i = PREF2/(rij*rij*rij)*mag_mu*mag_mu*(3 * (mixn*nx + miyn*ny) * (mjxn*ny - mjyn*nx) + mixn*mjyn - miyn*mjxn);
							torque_j = PREF2/(rij*rij*rij)*mag_mu*mag_mu*(3 * (-mjxn*nx - mjyn*ny) * (-mixn*ny + miyn*nx) + mjxn*miyn - mjyn*mixn);

							// Force on the grain i
							disk[i].fx+=-Fx;
							disk[i].fy+=-Fy;
							disk[i].Mz+= torque_j;
							//disk[i].fnorm+=sqrt((Fn*nx+ftx)*(Fn*nx+ftx) + (Fn*ny+fty)*(Fn*ny+fty));

							// Force on the grain j
							//disk[j].fx+=Fx;
							//disk[j].fy+=Fy;
							//disk[j].Mz+= torque_i;

							//print torques in console
							//if(bcltps%PASTPS_EPS==0){
							//		printf("t = %f, i=%d,j=%d\tr = (%.3f,%.3f)\tmin = (%.3f,%.3f)\tmjn = (%.3f,%.3f)\nscali=%.3f\tscalj=%.3f\tveci=%.3f\tcevj=%.3f\ntorque_i = %.3f, torque_j = %.3f\n", bcltps*DT, i, j, nx, ny, mixn,miyn,mjxn,mjyn, (mixn*nx + miyn*ny), (-mjxn*nx - mjyn*ny),(mjxn*ny - mjyn*nx), (-mixn*ny + miyn*nx) ,torque_i, torque_j);
							//}

							//print forces in console
							//if(bcltps%PASTPS_EPS==0){
							//		printf("t = %f, Ozi=%f,Ozj=%f\tr = (%.3f,%.3f)\tmin = (%.3f,%.3f)\tmjn = (%.3f,%.3f)\nFX = %f : term1=%.3f\tterm2=%.3f\tterm3=%.3f\tterm4=%.3f\nFY = %f : term1=%.3f\tterm2=%.3f\tterm3=%.3f\tterm4=%.3f\n", bcltps*DT, disk[i].Oz, disk[j].Oz, nx, ny, mixn,miyn,mjxn,mjyn, Fx, PREF1/(rij*rij*rij*rij)*(ny*mix - nx*miy) * mjy, PREF1/(rij*rij*rij*rij)*(ny*mjx - nx*mjy) * miy,- PREF1/(rij*rij*rij*rij)*2*nx*(mix*mjx + miy*mjy) , PREF1/(rij*rij*rij*rij)*5*nx*(ny*mix - nx*miy)*(ny*mjx - nx*mjy) ,Fy, PREF1/(rij*rij*rij*rij)*(nx*miy - ny*mix) * mjx, PREF1/(rij*rij*rij*rij)*(nx*mjy - ny*mjx) * mix, -PREF1/(rij*rij*rij*rij)*2*ny*(mix*mjx + miy*mjy),PREF1/(rij*rij*rij*rij)*5*ny*(ny*mix - nx*miy)*(ny*mjx - nx*mjy));
							//}


							//disk[j].Mz+= disk[j].Ray * (nx*fty-ny*ftx);
							//disk[j].fnorm+=sqrt((Fn*nx+ftx)*(Fn*nx+ftx) + (Fn*ny+fty)*(Fn*ny+fty));
						}
					}

				j=list_link[j];
				
				}
			}
		}

	}
}