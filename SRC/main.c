/* Copyright (c) China University of Petroleum, 2023.*/
/* All rights reserved.								*/

/* SUDEMOD2: $vision: 1.0 $ ; $Date: 2023/02/18  $	*/

#include "par.h"
#include "su.h"
#include <stdio.h>

#include "DemRun.h"
#include "Element.h"

/*********************** self documentation **********************/
char *sdoc[] = {
	" 									",
	" SUDEMOD2 - Discrete-Element MODeling (2nd order) for acoustic wave",
	"    equation with PML absorbing boundary conditions.                   ",
	"   ",
	" SUDEMOD2 [required parameters] [optional parameters]",
	" Required Parameters:							",
	"   ",
	" Optional Parameters:							",
	"   ",
	NULL};
/*
 * Authors: pengke, 2023, China University of Petroleum
 *
 * 
 * 
 */
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	int i, j;

	/* hook up getpar to handle the parameters */
	initargs(argc, argv);
	requestdoc(0);

	/* get required parameters */
	LoadParameters();

	/* Initialize data */
	GetSize();
	InitializeVar();

	/* Input and output file information */

	GetMesh();

	SetElementArguments(iModel, m_oBulkMaterial);
	SetContactArguments(iModel, m_oBulkMaterial);

	Preprocess();

	Calculation();

	/* free space before returning */
	FreeVar();

	printf("\n");
	printf("================SUDEMOD2 Program Running End================\n");
	return(CWP_Exit());
}

