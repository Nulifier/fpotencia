/*!
 * \file
 * \author Eric MSP Veith <eveith+fpotencia@veith-m.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef TYPES_H
#define TYPES_H


#include <complex>


namespace fPotencia {


    /*!
     * \brief The ValueType enum denotes whether a value is given in p.u. or
     *  as absolute value
     */
    enum ValueType
    {
        pu,         //!< Relative value in the per-unit system
        absolute    //!< Absolute value
    };


    //! \brief The type of the complex number used in fPotencia
    typedef std::complex<double> cx_double;
}


#endif // TYPES_H
