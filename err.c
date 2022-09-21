#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "def.h"


void  error
//==========================================================
//
//  PRINT ERROR MESSAGE AND EXIT THE PROGRAM
//
//
(
    int       error_type
)
//----------------------------------------------------------
{
    char      message[256];

    if (error_type == 0) sprintf(message,"NO POST PROCESS\n");
    if (error_type == 1) sprintf(message,"MALLOC ERROR\n");
    if (error_type == 2) sprintf(message,"FILE OPEN ERROR\n");
    if (error_type == 3) sprintf(message,"ERROR OF NORMAL UNIT VECTOR\n");
    if (error_type == 4) sprintf(message,"ERROR OF SORTING OF NODE LINK\n");
    if (error_type == 5) sprintf(message,"ERROR OF CALCULATION FOR ADHERENT SPRING\n");
    if (error_type == 6) sprintf(message,"ERROR OF DATA ARRAY\n");
    if (error_type == 7) sprintf(message,"ERROR OF DATA:DATA IS NEGATIVE\n");

    printf("%s\n",message);
    exit(1);

    return;
}

