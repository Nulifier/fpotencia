/*!
 * \file
 * \author Santiago Peñate Vera
 * \author Eric MSP Veith <eveith+fpotencia@veith-m.de>
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LOAD_H
#define	LOAD_H


#include <string>


#include "Bus.h"
#include "fpotencia_libs.h"


namespace fPotencia {


    /*!
     * \brief The Load class signifies a load that draws real and reactive
     *  power from the grid.
     */
    class Load {
    public:


        //! \brief Creates a new load object and initializes all values to 0
        explicit Load();


        //! \brief Returns the power as complex number
        cx_double power() const { return power_; }


        //! \brief Sets the load's power values
        Load& power(cx_double const& power)
        {
            power_ = power;
            return *this;
        }


        //! \brief Returns the real power drawn by the load.
        double realPower() const { return power_.real(); }


        //! \brief Sets the real power drawn by the load
        Load& realPower(double power)
        {
            power_.real(power);
            return *this;
        }


        //! \brief Returns the reactive power drawn by the load
        double reacivePower() const { return power_.imag(); }


        //! \brief Sets the reactive power drawn by the load
        Load& reactivePower(double power)
        {
            power_.imag(power);
            return *this;
        }


        //! \brief Connects the load to a bus
        Load& bus(Bus const& bus)
        {
            busIndex_ = bus.index;
            return *this;
        }


        size_t busIndex() const
        {
            assert(busIndex_ >= 0);
            return busIndex_;
        }


        //! \brief Retrieves the node's name
        std::string const& name() const
        {
            return name_;
        }


        //! \brief Sets a describing name of the node (for pretty-printing)
        Load& name(std::string const& name)
        {
            name_ = name;
            return *this;
        }


    private:


        //! \brief The name, for readability
        std::string name_;


        /*!
         * \brief The power this load draws.
         *
         * The real part of the complex number designates the real power
         * consumed, the imaginary part the reactive power that is consumed.
         */
        cx_double power_;


        //! \brief Index of the bus we're connected to.
        ptrdiff_t busIndex_;
    };
}

#endif	/* LOAD_H */
