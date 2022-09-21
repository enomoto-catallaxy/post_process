#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "def.h"

void  gpu_init (int argc, char **argv);

int  main(int argc, char *argv[])
{
   process   prc;
   domain    cdo;
   cell      cel;
   lattice   ltc;
   Wall      wall;

   int       flu_skip = 50,
	     skip  = 1,
             icnt  = 0,
	     max_icnt = 20000,
	     end_step = 0,
             count_flu;
   char      filename[256];

   // Set file point
   char      fpoint[256] = "./../";
 	 char 	 fname[256]="test";
   

//==============================================================================================
// Post process number  
// 1 :BIN to DAT(cell)
// 2 :Cell velocity
//==============================================================================================

   // Read configuration file
   sprintf(filename,"%s%s/result/config.dat",fpoint,fname);
   read_parameter(&prc,&cdo,&ltc,&cel,&wall,filename);

   while (icnt < max_icnt && end_step < 1) {
     // --- 1 --- Converting bin data to dat data(cell)
     end_step = BinToDat_Cell(icnt,&cdo,&ltc,&cel,&wall,fpoint,fname);
 
     // --- 2 --- 
     // end_step = Velocity(icnt,&prc,&cdo,&ltc,&cel,&wall,fpoint,fname);

     icnt += skip;
   }

   return 0;
}

