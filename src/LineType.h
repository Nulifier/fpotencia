/*!
 * \file
 * \author Eric MSP Veith <eveith+fpotencia@veith-m.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef FPOTENCIA_LINETYPE_H
#define FPOTENCIA_LINETYPE_H


#include "Types.h"


namespace fPotencia {


    /*!
     * \brief The LineType class defines types of lines given their series
     *  parameters R, X, and B.
     */
    class LineType
    {
    public:

        //! \brief Initializes all values to 0 (execpt b = 1e-9)
        explicit LineType():
                resistance_(0.0),
                reactance_(0.0),
                susceptance_(1e-9),
                valueType_(fPotencia::pu)
        {
        }


        virtual ~LineType() noexcept {}


        //! \brief Returns the line type's series resistance
        double resistance() const { return resistance_; }


        //! \brief Sets the line type's resistance
        LineType& resistance(double r)
        {
            resistance_ = r;
            return *this;
        }


        //! \brief Returns the line type's series reactance
        double reactance() const { return reactance_; }


        //! \brief Sets the line type's reactance
        LineType& reactance(double x)
        {
            reactance_ = x;
            return *this;
        }


        //! \brief Returns the line type's susceptance
        double susceptance() const { return susceptance_; }


        //! \brief Sets the line type's series susceptance (min 1e-9)
        LineType& susceptance(double b)
        {
            susceptance_ = (1.0 == b + 1.0 ? 1e-9 : b);
            return *this;
        }


        //! \brief Returns whether the values are absolute or in p.u.
        ValueType valueType() const { return valueType_; }


        //! \brief Sets whether the values are absolute or in p.u.
        LineType& valueType(ValueType type)
        {
            valueType_ = type;
            return *this;
        }


        //! \brief Returns the line type's complex impedance
        cx_double impedance() const
        {
            return cx_double(resistance(), reactance());
        }


        //! \brief Returns the line type's series shunt admittance
        cx_double shuntAdmittance() const
        {
            return cx_double(0, susceptance());
        }


    private:


        double resistance_;
        double reactance_;
        double susceptance_;
        ValueType valueType_;
    };
} // namespace fPotencia

#endif // FPOTENCIA_LINETYPE_H
