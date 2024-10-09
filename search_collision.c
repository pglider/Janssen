#include "decla.h"

void search_collision(void) {
	int IXcell, JYcell, ixTemp, jyTemp, ixTempbis;
	int contactcourant, bouclecontact;

	for(bouclecontact=0;bouclecontact<NMAXcontacts;bouclecontact++)
		  { 
		    disk[i].contactpr[bouclecontact]=disk[i].contact[bouclecontact];
		    disk[i].contact[bouclecontact]=-1;
		    disk[i].utijxpr[bouclecontact]=disk[i].utijx[bouclecontact];
		    disk[i].utijypr[bouclecontact]=disk[i].utijy[bouclecontact];
		    disk[i].utijx[bouclecontact]=0.;
		    disk[i].utijy[bouclecontact]=0.;
		  }
			
		// identify the cell part. i belongs to
		IXcell=(int) (disk[i].x/(grid*R));
		JYcell=(int) (disk[i].y/(grid*R));

		if( JYcell==NCELLY-1)
		printf("pb!: NCellY trop petit\t %d\t%d\n", JYcell, NCELLY);


		// look for possible collisions with part. in the neighbor cells
		for(ixTemp=IXcell-1;ixTemp<=IXcell+1 ; ixTemp++){
			for(jyTemp=JYcell-1;jyTemp<=JYcell+1 ; jyTemp++){
			   	ixTempbis=CLPcell(ixTemp,NCELLX);
				if(ixTempbis>=0 && jyTemp>=0){

				j=HoC[ixTempbis][jyTemp];// j = Head of Chain
				while(j!=0){// j=0 means end of chain
					if(i!=j)
					{

						// collision ?
						if( (CLPdist((disk[i].x -disk[j].x),WIDTH)*CLPdist((disk[i].x -disk[j].x ),WIDTH) + (disk[i].y - disk[j].y )*(disk[i].y -disk[j].y)) < ((disk[i].Ray + disk[j].Ray)*(disk[i].Ray + disk[j].Ray)))
							{  

							// Looking if the contact is already existing
							for(bouclecontact=0;bouclecontact<NMAXcontacts;bouclecontact++)
								{
								if(disk[i].contactpr[bouclecontact]==j)//existing contact
									{
									disk[i].utijx[disk[i].Nb_Contact]=disk[i].utijxpr[bouclecontact];
									disk[i].utijy[disk[i].Nb_Contact]=disk[i].utijypr[bouclecontact];
									break;
									}
								}
						
							//New contact
							contactcourant = disk[i].Nb_Contact;

							if(contactcourant>=NMAXcontacts) printf("\n ERREUR !!\n");
							disk[i].contact[contactcourant]=j;
							disk[i].Nb_Contact++;
						
							if(disk[i].Nb_Contact>=NMAXcontacts) printf("\n ERREUR i= %d t = %d %d \n",i,bcltps,disk[i].Nb_Contact);

							force_collision(contactcourant);
						}
					}

				j=list_link[j];
				
				}
			}

		}
	}
}