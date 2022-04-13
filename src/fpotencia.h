/* 
 * File:   fpotencia.h
 * Author: spv86_000
 *
 * Created on 31 de enero de 2015, 15:22
 */
#pragma once

/*Circuits*/
#include "Circuit.h"

/*Devices*/
#include "Bus.h"
#include "Load.h"
#include "Line.h"
#include "Transformer.h"
#include "Generator.h"
#include "Shunt.h"

/*Device 'designers'*/
#include "LineConstructors.h"
#include "TransformerConstructors.h"

/*Circuit solvers*/
//#include "Solver_ZBusGS.h"
//#include "Solver_NRrect.h"
//#include "Solver_Jacobi.h"
#include "Solver_NRpolar.h"
#include "Solver_NRcurrent.h"
#include "CommonFunctions.h"
