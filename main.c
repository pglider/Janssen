#include "decla.h"
#include "init.c"
#include "cell.c"
#include "search_collision.c"
#include "force_mag.c"
#include "createps.c"
#include "force_collision.c"
#include "force_ext.c"
#include "integre_Verlet.c"
#include "save_data.c"
#include "force_walls.c"

int main(){


    printf("No arguments given, using default folder OUTPUT and theta = 25 deg\n");
    //create save directory
    sprintf(save_folder,"OUTPUT/");
    mkdir(save_folder,0700);
    printf("Folder %s created\n",save_folder);
    //Clear eps in folder
    system("del OUTPUT\\*.eps");
    system("del OUTPUT\\*.jpg");
    theta = 0;

    
    //Print some parameters
    printf("N_PART = %d\n", N_PART);
    //initialize files
    printf("Initializing files...\n");
    sprintf(name_force,"force.dat");
    files_initialized = init_files();
    if(files_initialized){
        printf("Files initialized\n");
    }
    else{
        printf("Error initializing files\n");
        return 1;
    }
    // Load walls
    printf("Initializing walls...\n");
    load_walls();
    walls[0].x1 = 0;
    walls[0].y1 = Y0_PISTON;
    walls[0].x2 = WIDTH;
    walls[0].y2 = Y0_PISTON;
    walls[1].x1 = WIDTH;
    walls[1].y1 = 0;
    walls[1].x2 = WIDTH;
    walls[1].y2 = HEIGHT;
    walls[2].x1 = 0;
    walls[2].y1 = 0;
    walls[2].x2 = 0;
    walls[2].y2 = HEIGHT;
    eps_count=1;
    //initialize grains
    y_piston = Y0_PISTON;
    printf("Initializing grains...\n");
    init_grains();
    printf("Grains initialized\n NB_ITERATION = %d\n", NB_ITERATION);

    cell();
    
    for(bcltps=1;bcltps<=NB_ITERATION;bcltps++) {
        //Reset forces
        fy_piston = 0;
        fy_left = 0;
        fy_right = 0;
        
        //Update piston position
        if(bcltps==SEDIMENTATION_TIME){
            for(i=1;i<=N_PART;i++){
                disk[i].Oz=0;
                disk[i].Ozold=0;
                disk[i].dOz=0.;
            }
            printf("--------------------------=============PISTON STARTING=========-----------------\n");
        }
        if(bcltps>SEDIMENTATION_TIME){
            y_piston = Y0_PISTON - V_PISTON * (bcltps-SEDIMENTATION_TIME) * DT;
        }
        // Update the position of the first wall element with y_piston
        walls[0].y1 = y_piston;
        walls[0].y2 = y_piston;
        //Update magnetization parameter
        //Magnetization is set to zero during the sedimentation phase
        //Then it is ramped up linearly up to MAG_MU during the magnetization ramp time
        if(bcltps<SEDIMENTATION_TIME){
            mag_mu = 0.;
        }
        else{
            if(bcltps<SEDIMENTATION_TIME+MAGNETIZATION_RAMP_TIME){
                mag_mu = (bcltps-SEDIMENTATION_TIME)*MAG_MU/MAGNETIZATION_RAMP_TIME;
            }
            else{
                mag_mu = MAG_MU;
            }
        }
        //Initialize grid
        cell();

        //Clear force.dat file if needed
		if((DRAWFORCES)&&(bcltps%PASTPS_EPS==0)){
            fclose(fforce);
			fforce=fopen(name_force,"w");
            //return error if file cannot be opened
            if (fforce == NULL) {
                fprintf(stderr, "Error: Could not open file %s for writing, errno = %d\n", name_force, errno);
                return 1;
            }
		}


        // Loop over all grains
        for(i=1;i<=N_PART;i++){
            disk[i].highlight=0;
            //First the gravity field is taken into account
            force_ext();
            //Then collisions are detected and forces are computed
            disk[i].Nb_Contact=0;
            search_collision();
            if(disk[i].fixed==0){force_walls();}
            //Then the magnetic forces are computed
            if(disk[i].ismag==1){
                force_mag(mag_mu);
            }

        }

        //Knowing the forces, the equations of motion are integrated using Verlet's algorithm
        for(i=1;i<=N_PART;i++){

            if(disk[i].fixed==0){
                //Integrate equations of motion
                integre_Verlet();
            }
            else{
                //Do not integrate equations of motion for the inactive grains
                integre_VerletStop();
            }
        }
        

        //Save data every 100 timesteps
        if(bcltps%(ENERGY_FREQ)==1){

            //Reset variables
            Ec_tot = 0.;//Kinetic energy

            for(i=1;i<=N_PART;i++) {
                //compute kinetic energy
                Ec=0.5*disk[i].mass*disk[i].dx*disk[i].dx + 0.5*disk[i].mass*disk[i].dy*disk[i].dy  + 0.5*disk[i].II*disk[i].dOz*disk[i].dOz;
                if(Ec != 0){
                    Ec_tot+=Ec;//Add it to the total
                }
            }

            //Rescale Ec_tot
            Ec_tot = Ec_tot/(double)N_PART;

            //Write energy to files
            fprintf(fp_energy, "%d\t%e\n", bcltps, Ec_tot);
        }

        if((bcltps%PASTPS_EPS==0)&&(bcltps<PASTPS_EPS_END)&&(bcltps>SEDIMENTATION_TIME)){
            createps(save_folder,eps_count);
            eps_count = eps_count + 1;
            printf("Image number %d/%d created, Kinetic energy = %e, fy_piston = %8f, fy_left = %8f, fy_right = %8f\n", eps_count, (PASTPS_EPS_END-SEDIMENTATION_TIME)/PASTPS_EPS, Ec_tot, fy_piston/(10*MASS), fy_left/(10*MASS), fy_right/(10*MASS));
        }

        //In case of error such as high interpenetration, write an image to see where the error occured
        if(err_count!=last_err_count){
            err=1;
            createps(save_folder, bcltps/PASTPS_EPS);
        }

        //End of simulation
        if(err_count>10){
            //save final state
            printf("\n\n[WARNING] Exiting after too many errors, bcltps = %d\n",bcltps);
            err=1;
            createps(save_folder, bcltps/PASTPS_EPS);
            //Write final pos file and image
            bcltps = bcltps + POS_FREQ;//Make sure a new pos file is created
            save_data();
            break;
        }

        //Save positions
        if(bcltps%POS_FREQ==0){
            save_data();
        }

        //Update and reset error counter
        last_err_count=err_count;
        err=0;

        //Print progression every 5%
        if(bcltps%(NB_ITERATION/20)==0){
            printf("progress : %d/100 --- log10(Ec) = %.2f\n",5*bcltps/(NB_ITERATION/20),log10(Ec_tot));
        }

    }

    //End of simulation
    //save final image
    printf("End of simulation at iteration %d\n",bcltps-1);
    createps(save_folder, bcltps/PASTPS_EPS);

    //Close files and end program
    fclose(fp_energy);
    return 0;
}
