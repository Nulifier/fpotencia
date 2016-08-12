/*!
 * \file
 * \author Santiago Peñate Vera
 * \author Eric MSP Veith <eveith+fpotencia@veith-m.de>
 *
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef LINE_H
#define	LINE_H


#include <cstddef>

#include "Types.h"
#include "LineType.h"
#include "Solution.h"
#include "fpotencia_libs.h"


using namespace std;

namespace fPotencia {
    class Circuit;


    /*!
     * \brief The Line class represents an actual power line
     */
    class Line {
    public:
        explicit Line():
                length_(0),
                bus1Id_(-1),
                bus2Id_(-1),
                impedance_(0, 0),
                shuntAdmittance_(0, 0)
        {
        }


        virtual ~Line() noexcept {};


        /*!
         * \brief Returns the type of the line
         *
         * \sa LineType
         */
        LineType const& lineType() const { return lineType_; }


        /*!
         * \brief Sets the new line type
         *
         * Setting the line type forces the recalculation of the line's
         * values.
         *
         * \param[in] type The new line type
         *
         * \return `*this`
         */
        Line& lineType(LineType const& type)
        {
            lineType_ = type;
            calculateLineParameters();
            calculateAdmittance();
            return *this;
        }


        //! \brief Returns the line's length
        double length() const { return length_; }


        /*!
         * \brief Sets the line's length
         *
         * Setting the length of the line forces the recalculation of the
         * line's parameters.
         *
         * \param[in] length The new line length
         *
         * \return `*this`
         */
        Line& length(double length)
        {
            length_ = length;
            calculateLineParameters();
            calculateAdmittance();
            return *this;
        }


        //! \brief Returns the complex impedance of the line
        cx_double impedance() const { return impedance_; }


        //! \brief Sets the line's impedance and recalculates the admittance
        Line& impedance(cx_double const& impedance)
        {
            impedance_ = impedance;
            calculateAdmittance();
            return *this;
        }


        //! \brief Returns the complex impedance of the line
        cx_double shuntAdmittance() const { return shuntAdmittance_; }


        /*!
         * \brief Sets the line's shunt admittance and recalculates
         *  the admittance
         */
        Line& shuntAdmittance(cx_double const& shuntAdmittance)
        {
            shuntAdmittance_ = shuntAdmittance;
            calculateAdmittance();
            return *this;
        }


        //! \brief Returns the type of values: Absolute or in p.u.
        ValueType valueType() const { return valueType_; }


        //! \brief Sets the value type: absolute or in p.u.
        Line& valueType(ValueType valueType)
        {
            valueType_ = valueType;
            return *this;
        }


        //! \brief Retrieves the base value for the impedance (if p.u.)
        double impedanceBase() const { return impedanceBase_; }


        //! \brief Returns the IDs of the two buses the line connects
        std::pair<size_t, size_t> buses() const
        {
            assert(bus1Id_ >= 0);
            assert(bus2Id_ >= 0);
            return std::make_pair(bus1Id_, bus2Id_);
        }


        //! \brief Connects two buses with the line
        Line& buses(size_t bus1Id, size_t bus2Id)
        {
            bus1Id_ = bus1Id;
            bus2Id_ = bus2Id;
            return *this;
        }


        /*!
         * \brief Places the line's admittance values in the system-wide
         *  admittance matrix
         *
         * \param[inout] systemAdmittance The admittance matrix of the whole
         *  circuit
         */
        void placeAdmittance(sp_cx_mat& systemAdmittance) const;


        /*!
         * \brief Sets the power flowing between the buses on this line
         *
         * \param[in] solution A circuit solution as calculated by a solver
         */
        void setPower(cx_solution const& solution);


        /*!
         * \brief Returns the power that flows between the buses.
         *
         * \return A `std::pair` containing (1) the power flowing from bus 1
         *  to bus 2, and (2) the power flowing from bus 2 to bus 1
         */
        std::pair<cx_double, cx_double> power() const
        {
            return std::make_pair(powerFlowBus1ToBus2_, powerFlowBus2ToBus1_);
        }


        //! \brief Returns the power loss
        cx_double powerLoss() const
        {
            return powerLoss_;
        }


        /*!
         * \brief Returns the current that flows between the buses.
         *
         * \return A `std::pair` containing (1) the current flowing from bus 1
         *  to bus 2, and (2) the current flowing from bus 2 to bus 1
         */
        std::pair<cx_double, cx_double> current() const
        {
            return std::make_pair(
                    currentFlowBus1ToBus2_,
                    currentFlowBus2ToBus1_);
        }


    private:


        friend class Circuit;


        double length_;
        LineType lineType_;
        std::ptrdiff_t bus1Id_;
        std::ptrdiff_t bus2Id_;

        ValueType valueType_;
        cx_double impedance_;
        double impedanceBase_;
        cx_double shuntAdmittance_;


        //! \brief Calculated 2x2 admittance matrix
        cx_mat admittance_;

        cx_double powerLoss_;
        cx_double powerFlowBus1ToBus2_;
        cx_double powerFlowBus2ToBus1_;
        cx_double currentFlowBus1ToBus2_;
        cx_double currentFlowBus2ToBus1_;

        //! \brief Calculates the line's parameters
        void calculateLineParameters();

        //! \brief Re-calculates the line's admittance matrix
        void calculateAdmittance();
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
