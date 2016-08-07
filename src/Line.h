/*!
 * \file
 * \author Santiago Peñate Vera
 *
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef LINE_H
#define	LINE_H


#include "LineType.h"
#include "Solution.h"
#include "fpotencia_libs.h"


using namespace std;

namespace fPotencia {


    /*!
     * \brief The Line class represents an actual power line
     */
    class Line {
    public:
        Line(string name, int connection_bus1, int connection_bus2, LineType line_type, double line_lenght);

        void SetType(LineType line_type);

        virtual ~Line();

        //properties
        string Name;

        int bus1 = 0;

        int bus2 = 0;

        double lenght = 0;

        bool values_in_per_unit;

        void get_element_Y(int n, sp_cx_mat &Yret);

        void calculate_current(cx_solution sol);

        void print();


        /*************************************************************************
         * Calculated variables: Results
         *************************************************************************/

        cx_double current_bus1_to_bus2;

        cx_double current_bus2_to_bus1;

        cx_double power_bus1_to_bus2;

        cx_double power_bus2_to_bus1;

        cx_double power_losses;

        double Zbase;

    private:

        cx_double impedance;

        cx_double shunt_admittance;

        /*************************************************************************
         * Calculated variables: Results
         *************************************************************************/
        cx_mat Y_element; //calculated element admittance matrix (2x2)
    };

    /*******************************************************************************
     *Line type for 3-phase lines
     ******************************************************************************/
    class LineType3 {
    public:
        LineType3(string name, cx_mat3 Z_abc, cx_mat3 Y_abc);

        virtual~LineType3();


        //properties
        string Name;

        cx_mat3 Zabc;

        cx_mat3 Yabc;

    private:

    };

    /*
     *Line for 3-phase
     */
    class Line3 {
    public:
        Line3(string name, int connection_bus1, int connection_bus2, LineType3 line_type, double line_lenght);

        virtual~Line3();

        void SetType(LineType3 &line_type);

        //properties
        string Name;

        int bus1 = 0;

        int bus2 = 0;

        double lenght = 0;


    private:

        cx_mat3 A;
        cx_mat3 B;

        cx_mat3 a;
        cx_mat3 b;
        cx_mat3 c;
        cx_mat3 d;
    };
}

#endif	/* LINE_H */
