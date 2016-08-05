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


#ifndef GENERATOR_H
#define	GENERATOR_H


#include <string>

#include "Bus.h"
#include "fpotencia_libs.h"


namespace fPotencia {
    class Bus;


    /*!
     * \brief The Generator class represents a generator in a power grid.
     */
    class Generator
    {
    public:


        //! \brief Unit type of the voltage value
        enum VoltageType
        {
            pu,     //!< Per Unit
            volts   //!< Volts absolute
        };


        //! \brief Creates a new, unitialized, unconnected Generator object
        explicit Generator();


        //! \brief Returns the generator name
        std::string const& name() const { return name_; }


        //! \brief Sets the name of the generator, for prettyprinting.
        Generator& name(std::string const& name)
        {
            name_ = name;
            return *this;
        }


        //! \brief Connects the generator to a bus.
        Generator& bus(Bus const& bus)
        {
            busIndex_ = bus.index;
            return *this;
        }


        /*!
         * \brief Returns the Bus the generator is connected to.
         *
         * If the generator has not been connected to a bus yet, this will
         * return a `nullptr`!
         *
         * \return The Bus object, may be `nullptr`
         */
        size_t busIndex() const
        {
            assert(busIndex_ > 0);
            return busIndex_;
        }


        //! \brief The voltage the generator supplies to the grid
        double voltage() const { return voltage_; }


        /*!
         * \brief Returns the unit type of the voltage: P.u. or absolute volts
         *
         * \sa VoltageType
         */
        VoltageType voltageType() const { return voltageType_; }


        /*!
         * \brief Sets the voltage supplied by the generator
         *
         * \param voltage The voltage value
         *
         * \param voltageType The type of the voltage value: P.u. or absolte
         *  volts
         *
         * \return `*this`
         */
        Generator& voltage(double voltage, VoltageType voltageType)
        {
            voltage_ = voltage;
            voltageType_ = voltageType;
            return *this;
        }


        //! \brief Returns the real power of the generator
        double realPower() const { return power_.real(); }


        //! \brief Sets the amount of real power supplied by the generator
        Generator& realPower(double p)
        {
            power_.real(p);
            return *this;
        }


        /*!
         * \brief Returns the reactive power suppled by the generator
         *
         * The reactive power supplied by the generator is only available
         * when the power system load flow analysis has been conducted.
         *
         * \sa Solver
         */
        double reactivePower() const { return power_.imag(); }


        //! \brief Returns the complex power value of this generator
        cx_double power() const { return power_; }


        //! \brief Sets the VAr limits for the generator
        Generator& reactivePowerLimits(double min, double max)
        {
            minQ_ = min;
            maxQ_ = max;
            return *this;
        }


        //! \brief Returns the minimum and maximum reactive power limit
        std::pair<double, double> reactivePowerLimits() const
        {
            return std::make_pair(minQ_, maxQ_);
        }


    private:



        //! \brief Optional name of the generator, for pretty-printing.
        std::string name_;


        //! \brief The index of the Bus we're connected to
        ptrdiff_t busIndex_;


        //! \brief Voltage supplied by the generator
        double voltage_;


        //! \brief Is ::voltage_ in p.u. or volts absolute?
        VoltageType voltageType_;


        /*!
         * \brief The power supplied by the generator
         *
         * The real part of the complex number signifies the amount of
         * real power provided; the imaginary part the amount of reactive
         * power fed in by the generator.
         */
        cx_double power_;


        //! \brief The minimum VAr output of the generator
        double minQ_;


        //! \brief Upper limit of reactive power output of the generator
        double maxQ_;
    };
}

#endif	/* GENERATOR_H */
