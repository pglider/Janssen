#include "decla.h"

void createps(char * save_folder, int ieps){
	char name[255],fileps[255]="",str[50];
	double xref,xref2,yref,xi,yi,xj,yj,fn;
	int forcetype;

	if(err==1){
        sprintf(name,"err%.5d.eps",ieps);
	}
	else{
		sprintf(name,"p%.5d.eps",ieps);
	}

	if(bcltps>=NB_ITERATION){
		sprintf(name,"final.eps");
	}

	//Open image file
	strcat(fileps,save_folder);
	strcat(fileps,name);
	//Try to open the file, or else produce an error
	pos_eps = fopen(fileps, "w");
	if (pos_eps == NULL) {
		fprintf(stderr, "Error: Could not open file %s for writing\n", fileps);
		return;
	}

	//printf("Calling createps with name %s, fileps %s, and pos_eps %p\n", name, fileps, pos_eps);

	fprintf(pos_eps,"%%!PS-Adobe-3.0 EPSF-3.0\n");

    //Bounding box
	fprintf(pos_eps,"%%");
	//sprintf(str,"%%BoundingBox: %d %d %d %d\n",(int)(-2*R*SCALE),(int)(-2*R*SCALE),(int)((WIDTH+2*R)*SCALE),(int)(0.7*HEIGHT*SCALE));
	sprintf(str,"%%BoundingBox: %d %d %d %d\n",(int)(-2*R*SCALE),(int)(0.25*HEIGHT*SCALE),(int)((WIDTH+2*R)*SCALE),(int)(0.5*HEIGHT*SCALE));
	//Zoom on the first wall
	//sprintf(str,"%%BoundingBox: %d %d %d %d\n",(int)(2*R*SCALE),(int)(0.031*SCALE),(int)((25*R)*SCALE),(int)(0.037*SCALE));
	
	fprintf(pos_eps,"%s",str);
	fprintf(pos_eps,"%%");
	fprintf(pos_eps,"%%");
	fprintf(pos_eps,"%%EndComments \n");

	// Define shortcuts
	fprintf(pos_eps,"/M{moveto}def\n");//move to
	fprintf(pos_eps,"/L{lineto stroke}def\n");//line to
	fprintf(pos_eps,"/A{newpath 0 360 arc closepath stroke}def\n");
	fprintf(pos_eps,"/B{arc closepath gsave grestore fill}def\n");// For magnetic particles
	fprintf(pos_eps,"/C{0.75 1 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");
	fprintf(pos_eps,"/D{0.3 1 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");
	fprintf(pos_eps,"/E{0.2 0.1 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore stroke}def\n");
	fprintf(pos_eps,"/F{0.2 0.1 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");
	fprintf(pos_eps,"/G{sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");
	fprintf(pos_eps,"/H{sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore stroke}def\n");

	fprintf(pos_eps,"100 setgray\n");

    //Background
    fprintf(pos_eps,"0 0 60000 0.55 0.76 0.22 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke\n");

	// Draw walls
	fprintf(pos_eps,"5 setlinewidth 0 1 0 sethsbcolor\n");
	for (int i = 0; i < 4; i++) {
		fprintf(pos_eps, "newpath\n");
		fprintf(pos_eps, "%.3f %.3f moveto\n", walls[i].x1 * SCALE, walls[i].y1 * SCALE);
		fprintf(pos_eps, "%.3f %.3f lineto\n", walls[i].x2 * SCALE, walls[i].y2 * SCALE);
		fprintf(pos_eps, "stroke\n");
	}

	fprintf(pos_eps,"1 setlinewidth\n");
    
	//Grains in the packing
	for(i=1;i<=N_PART;i++) {

		xref=disk[i].x;
		yref=disk[i].y;

		//Lower boundary grains
		if(disk[i].fixed){

			//Draw clones
			// if(disk[i].x>0.95*WIDTH){
			// 	xref2 = xref - WIDTH;
			// 	fprintf(pos_eps,"\n%.3f %.3f %.3f %d 1 %d H", xref2*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0,0);
			// }
			// if(disk[i].x<0.05*WIDTH){
			// 	xref2 = xref + WIDTH;
			// 	fprintf(pos_eps,"\n%.3f %.3f %.3f %d 1 %d H", xref2*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0,0);
			// }

			fprintf(pos_eps,"\n%.3f %.3f %.3f 0 0 0.65 G", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE);
			fprintf(pos_eps,"\n%.3f %.3f %.3f %d 1 %d H", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0,0);
		}
		else{
			xref=disk[i].x;
			yref=disk[i].y;
			
			//Draw clones
			// if(disk[i].x>0.9*WIDTH){
			// 	xref2 = xref - WIDTH;
			// 	fprintf(pos_eps,"\n%.3f %.3f %.3f %d 1 %d H", xref2*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0,0);
			// }
			// if(disk[i].x<0.1*WIDTH){
			// 	xref2 = xref + WIDTH;
			// 	fprintf(pos_eps,"\n%.3f %.3f %.3f %d 1 %d H", xref2*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0,0);
			// }

			//Color according to sum of norm of forces
			//fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f C", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0.12 + fmin(700*disk[i].fnorm,0.88));

			//Color according to instantaneous speed
			if(disk[i].highlight==0){
			    fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f 1 1 G", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, fmin(fmax(0.5 + disk[i].Oz*0.02,0.),1.));
			}
			else{			
			    fprintf(pos_eps,"\n%.3f %.3f %.3f 0 C", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE);
			}
			
			//Color according to number of contacts
			//fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f C", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0.1 + disk[i].Nb_Contact*0.4/6);

			//Color according to rotation angle
			//fprintf(pos_eps,"\n %f %f %f %f D", (xref)*SCALE, yref*SCALE, disk[i].Ray*SCALE, (3.2 - disk[i].Oz)/6.4);

			//Color according to average speed
			//fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f C", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0.12 + fmin(0.5*(VLENETH*disk[i].v_avg),0.88));

			//Color according to momentum (used for debugging divergence caused by rotational friction implementation)
			//fprintf(pos_eps,"\n %f %f %f %f C", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0.55 + fmin(log10(10+abs(disk[i].Mz))/10,0.45));

			//Circle around carried grains
	//            if(disk[i].carried==1){
	//                fprintf(pos_eps,"\n%.3f %.3f %.3f 0.15 E", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE);
	//            }

			fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f 1 %.3f H", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 0.55 + fmin(300*disk[i].fnorm,0.45),0.);
			fprintf(pos_eps,"\n%f %f M", xref*SCALE, yref*SCALE);
			fprintf(pos_eps,"\n%f %f L\n", SCALE * (xref + disk[i].Ray*cos(disk[i].Oz)),SCALE * (yref + disk[i].Ray*sin(disk[i].Oz)));


			//if particle is magnetic, draw half a circle black
			if(disk[i].ismag==1){
				fprintf(pos_eps,"\n%f %f %f %f %f B", xref*SCALE, yref*SCALE, disk[i].Ray*SCALE, 57.2958*disk[i].Oz-90, 57.2958*disk[i].Oz+90);
			}
			//(EPS_SCALE*part_pos_x(i)), (EPS_SCALE*part_pos_y(i)), (EPS_SCALE*particle_r(i)), 57.2958*particle_angle(i)-90, 57.2958*particle_angle(i)+90

			//Place a dot on moving grains
	//            if(disk[i].moving==1){
	//                fprintf(pos_eps,"\n%.3f %.3f %.3f 0.75 F", xref*SCALE, yref*SCALE, 0.5*disk[i].Ray*SCALE);
	//            }

			//Place a dot colored according to force
			//fprintf(pos_eps,"\n%.3f %.3f %.3f %.3f 1 %.3f G", xref*SCALE, yref*SCALE, 0.5*disk[i].Ray*SCALE, (3.2 - disk[i].Oz)/6.4,fmax(0,1 - fmin(4000*disk[i].fnorm,1)));
		}
	
	}

	fprintf(pos_eps,"\n0.5 setlinewidth");
    fprintf(pos_eps,"\n1 0 0.5 sethsbcolor");
	fprintf(pos_eps,"\n[5 5] 0 setdash\n");
    fprintf(pos_eps,"\n%f %.1f M", 0*SCALE, 0.);
    fprintf(pos_eps,"\n%f %f L", 0*SCALE, HEIGHT*SCALE);
    //right wall
    fprintf(pos_eps,"\n%f %.1f M", WIDTH*SCALE, .0);
    fprintf(pos_eps,"\n%f %f L", WIDTH*SCALE, HEIGHT*SCALE);

	fprintf(pos_eps,"\n[] 0 setdash\n");
	if(DRAWFORCES){
		fclose(fforce);
		fforce=fopen(name_force,"r");
		// Read forces from file
		char line[256];
		while (fgets(line, sizeof(line), fforce)) {
			sscanf(line, "%lf %lf %lf %lf %lf %d", &xi, &yi, &xj, &yj, &fn, &forcetype);
			switch (forcetype) {
				case 0://Normal forces between grains
					fprintf(pos_eps,"\n0.1 0.1 0.1 sethsbcolor");
					fprintf(pos_eps, "\n%f setlinewidth\n", fmax(fmin(fn * 50000, 20),1));
					fprintf(pos_eps, "%f %f M %f %f L\n", xi * SCALE, yi * SCALE, xj * SCALE, yj * SCALE);
					break;
				case 1://Tangential forces between grains
					fprintf(pos_eps,"\n1 0.6 0.8 sethsbcolor");
					fprintf(pos_eps, "\n%f setlinewidth\n", fmax(fmin(fn * 50000, 20),1));
					fprintf(pos_eps, "%f %f M %f %f L\n", xi * SCALE, yi * SCALE, xj * SCALE, yj * SCALE);
					break;
				case 2://Normal forces grain-wall
					fprintf(pos_eps,"\n0.3 0.9 0.8 sethsbcolor");
					fprintf(pos_eps, "\n%f setlinewidth\n", fmax(fmin(fn * 50000, 20),1));
					fprintf(pos_eps, "%f %f M %f %f L\n", xi * SCALE, yi * SCALE, xj * SCALE, yj * SCALE);
					break;
				case 3://Tangential forces grain-wall
					fprintf(pos_eps,"\n0.2 0.3 1 sethsbcolor");
					fprintf(pos_eps, "\n%f setlinewidth\n", fmax(fmin(fn * 50000, 20),1));
					fprintf(pos_eps, "%f %f M %f %f L\n", xi * SCALE, yi * SCALE, xj * SCALE, yj * SCALE);
					break;
			}
		}
		fclose(fforce);
	}

	// Add text at the bottom of the image stating the value of fy
	// fprintf(pos_eps, "\n1 1 1 sethsbcolor\n");
	// fprintf(pos_eps, "%.0f %.0f moveto\n", 0.2 * WIDTH * SCALE, y_piston * SCALE);
	// fprintf(pos_eps, "(fy = %.2f) show\n", fy_piston);

	// //Add text at the top of the image stating the current time and angle
	// fprintf(pos_eps,"\n0 0 0 sethsbcolor\n");
	// fprintf(pos_eps, "%.0f %.0f moveto\n", 0.1*WIDTH, SCALE * 0.7 * HEIGHT + SCALE * 0.1 * WIDTH - SCALE/300);
	// fprintf(pos_eps,"(theta = %.1f, t = %.1f, log10(Ec) = %.2f) show\n",theta,bcltps*DT,log10(Ec_tot));

	// fprintf(pos_eps,"\nshowpage\n");

	fclose(pos_eps);

    //Optional : convert and delete to save disk space, I suggest not to do that
	//system("mogrify -format jpg p*eps\n");
    //system("rm p*.eps\n");
}
