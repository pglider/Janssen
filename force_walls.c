#include "decla.h"

void force_walls(void) {
    double x, y, rij;
    double nx, ny;
    double Fn;
    double Vn, Vs, Vsx, Vsy, vsx, vsy;
    double ftx, fty;
    double delta;
    double utx, uty, ut;
    int wall_index;
    //Loop over the walls
    for(wall_index=0;wall_index<4;wall_index++){
        //Compute the distance between the grain and the wall
        double x1 = walls[wall_index].x1;
        double y1 = walls[wall_index].y1;
        double x2 = walls[wall_index].x2;
        double y2 = walls[wall_index].y2;

        double dx = x2 - x1;
        double dy = y2 - y1;
        double t = ((disk[i].x - x1) * dx + (disk[i].y - y1) * dy) / (dx * dx + dy * dy);
        t = fmax(0, fmin(1, t));

        double closestX = x1 + t * dx;
        double closestY = y1 + t * dy;

        x = closestX - disk[i].x;
        y = closestY - disk[i].y;
        rij = sqrt(x * x + y * y);
        //Compute the normal vector
        nx = x / rij;
        ny = y / rij;
        //Compute the overlap
        delta = disk[i].Ray - rij;
        if(delta > 0){
            //printf("Overlap between grain %d (x=%f, y=%f) and wall %d (x1=%f, y1=%f, x2=%f, y2=%f) is %f\n", i, disk[i].x, disk[i].y, wall_index, x1, y1, x2, y2, delta);
            //disk[i].highlight=1;
            
            // Looking if the contact is already existing
            if(disk[i].contactwallpr[wall_index]==1){
                disk[i].utijx_walls[wall_index]=disk[i].utijxpr_walls[wall_index];
                disk[i].utijy_walls[wall_index]=disk[i].utijypr_walls[wall_index];
            }
            else{
                //New contact
                disk[i].contactwallpr[wall_index]=1;
                disk[i].utijx_walls[wall_index]=0.;
                disk[i].utijy_walls[wall_index]=0.;
            }

            
            // Compute the normal speed Vn=\vec{v}.\vec{n}, viscoelastic model
            Vn= - disk[i].dx*nx - disk[i].dy*ny;
            Fn=- ( KN * delta - GAMMA*Vn);
            
            //Compute the normal speed
            Vn = disk[i].dx * nx + disk[i].dy * ny;

            // Component of the sliding speed (relative speed - normal speed + Varignon formula), then norm
            Vsx= - disk[i].dx - Vn*nx + ny*disk[i].dOz*disk[i].Ray;
            Vsy= - disk[i].dy - Vn*ny - nx*disk[i].dOz*disk[i].Ray;
            Vs= sqrt(Vsx*Vsx+Vsy*Vsy);

            // Cundall's model u_t=\int{\vec{t}.\vec{dl}}
            utx=disk[i].utijx_walls[wall_index]+KT*Vsx*DT;
            uty=disk[i].utijy_walls[wall_index]+KT*Vsy*DT;
            ut=sqrt(utx*utx+uty*uty);
            if(fabs(ut/(MU_WALL*Fn))>1. && (utx*Vsx)>0 && (uty*Vsy)>0) {
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
                // The friction force is proportional to the normal force with a coefficient MU_WALL and is opposed to the sliding speed (Fn<0 if there is a collision, see above)
                ftx=- MU_WALL*Fn*vsx;
                fty=- MU_WALL*Fn*vsy;

                disk[i].utijx_walls[wall_index]=ftx;
                disk[i].utijy_walls[wall_index]=fty;
                disk[i].highlight=1;
            }
            else{
                disk[i].utijx_walls[wall_index]+=KT*Vsx*DT;
                disk[i].utijy_walls[wall_index]+=KT*Vsy*DT;
                ftx=disk[i].utijx_walls[wall_index];//+ GAMMA_PREFACTOR*GAMMA*reff/(.5*R)*Vsx;
                fty=disk[i].utijy_walls[wall_index];//+ GAMMA_PREFACTOR*GAMMA*reff/(.5*R)*Vsy;
            }

            //Debugging observed cohesion : Taking fn negative only
            //if(Fn>0){
            //    Fn=0;
            //}

            //Forces and moments
            disk[i].fx+=Fn*nx+ftx;
            disk[i].fy+=Fn*ny+fty;
            disk[i].Mz+= disk[i].Ray * (nx*fty-ny*ftx);

            switch(wall_index){
                case 0:
                    fy_piston+=Fn*ny;
                    break;
                case 1:
                    fy_left+=fty;
                    break;
                case 2:
                    fy_right+=fty;
                    break;
            }

            //Save the forces in a file
            if((DRAWFORCES)&&(bcltps%PASTPS_EPS==0)){
                fprintf(fforce,"%f\t%f\t%f\t%f\t%f\t%d\n",closestX,closestY,disk[i].x,disk[i].y,fabs(Fn),2);
                fprintf(fforce,"%f\t%f\t%f\t%f\t%f\t%d\n",closestX,closestY,closestX+R*0.1*ftx/MASS,closestY+R*0.1*fty/MASS,fabs(Fn),3);
            }

        }
        else{
            //No contact
            disk[i].contactwallpr[wall_index]=0;
        }
    }
}
