#include "decla.h"

int init_files() {
    //open energy.txt file to store kinetic energy data
    strcat(energyfilename,save_folder);
    strcat(energyfilename,"energy.txt");
    //Check energy.txt pointer is valid, this is to debug segmentation fault
    if ((fp_energy = fopen(energyfilename, "w")) == NULL) {
        fprintf(stderr, "Cannot open energy.txt");
        return 0;
    }
    return 1;
}

void init_grains(void) {
	double Ri, tcoll, Y0;
	int n,in;
	float Vo = 100*R;

	srand(SEED);

	// Display duration of a collision and restitution coefficient
    tcoll = 3.14/sqrt(KN/MASS - 0.25*(GAMMA*GAMMA)/(MASS*MASS));
    printf("Duration of a collision DT = %f = %.0f timesteps\n", tcoll,tcoll/DT);
    printf("Restitution e = %f\n\n",exp(-GAMMA/(2*MASS) * tcoll));

	Y0=y_piston+2*R;

	i=1;
    n=1;
	p0=i;
    // Rest of the particles
	while(i<=N_PART) {
    // We define a new layer in y
        // We fill the layer of thickness Rmax*n (Rmax=grid*R is the maximum radius a grain can have).
        // We add a random velocity (Vo initialized in line 8 of init), always to try to have a random stacking of grains
        for(in=0;in<(int)L;in++) {

            if((grid*R*(0.5 + in)>1*R)&&(grid*R*(0.5 + in)<WIDTH-1*R)){

                if(i%100==0){
                    printf("Particle %d out of %d\n", i, N_PART);
                }

                disk[i].x = grid*R*(0.5 + in) + 2*((double) rand()/RAND_MAX-0.5)*R*(DENSITE-2.4);
                disk[i].y = Y0 + grid*R*(n-1) + 2*((double) rand()/RAND_MAX-0.5)*R*(DENSITE-2.4);

                disk[i].xold = disk[i].x + Vo*DT*((double) rand()/RAND_MAX-0.5);;
                disk[i].yold = disk[i].y + Vo*DT*((double) rand()/RAND_MAX-0.5);
                // We define the radius of the grains: 1+0.2*(random number between -0.5 and 0.5), i.e. 0.9<r<1.1
                Ri= 1. + 2*POLYDISP*( (double) rand()/RAND_MAX-0.5 );
                // Initializations of the radii with the factors R, m and I defined in decla.h Pseudo-3D model (evolution in r² for the mass, in r⁴ for the moment of inertia).
                disk[i].Ray = R*Ri;
                disk[i].mass = MASS*Ri*Ri*Ri;
                disk[i].II = I*Ri*Ri*Ri*Ri;

                // Initialization of contacts: no contact with all particles j, the components (x, y, z) of the tangential displacement (ut) between i and j (ij) are zero.
                for(j=0;j<NMAXcontacts;j++) {
                    disk[i].contact[j]=-1;
                    disk[i].contactpr[j]=-1;
                    disk[i].utijx[j]=0.;
                    disk[i].utijy[j]=0.;
                }
                //Same for the walls
                for(j=0;j<MAXWALLS;j++) {
                    disk[i].utijx_walls[j]=0.;
                    disk[i].utijy_walls[j]=0.;
                }
                
                disk[i].dx=0.;
                disk[i].dy=0.;

                //Random initial orientation
                disk[i].Oz = 2*3.14*((double) rand()/RAND_MAX);
                disk[i].Ozold=disk[i].Oz;
                disk[i].dOz=0.;
                disk[i].fixed = 0;
                disk[i].highlight=0;

                i++;
            }

			if(i>N_PART){
                break;
            }
		}
		n++;
    }

    // Make NMAG random particles magnetic
    for(i=1;i<=NMAG;i++){
        j=1;
        while(disk[j].fixed || disk[j].ismag){
            j = 1+(int) ((N_PART) * ((double) rand()/RAND_MAX));
        }
        disk[j].ismag = 1;
    }

    //Print list of magnetic particles
    //printf("Magnetic particles : ");
    //for(i=1;i<=N_PART;i++){
    //    if(disk[i].ismag == 1){
    //        printf("%d ", i);
    //    }
    //}

}

void load_walls() {
    FILE *fp_walls;
    char line[256];
    int wall_index = 0;

    fp_walls = fopen("walls.txt", "r");
    if (fp_walls == NULL) {
        fprintf(stderr, "Cannot open walls.txt");
        return;
    }

    while (fgets(line, sizeof(line), fp_walls)) {
        printf("%s\n", line);
        if (sscanf(line, "%lf %lf %lf %lf", &walls[wall_index].x1, &walls[wall_index].y1, &walls[wall_index].x2, &walls[wall_index].y2) == 4) {
            printf("Wall %d: %f %f %f %f\n", wall_index, walls[wall_index].x1, walls[wall_index].y1, walls[wall_index].x2, walls[wall_index].y2);
            wall_index++;
        }
    }

    fclose(fp_walls);
}