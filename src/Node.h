/*!
 * \file
 * \author Eric MSP Veith <eveith+fpotencia@veith-m.de>
 *
 * Copyright (C) 2016 Eric MSP Veith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef FPOTENCIA_NODE_H
#define FPOTENCIA_NODE_H


#include <string>
#include <cstddef>

#include "Bus.h"


namespace fPotencia {

    /*!
     * \brief The Node class is the parent class of all nodes that can be
     *  connected to a bus such as loads and generators.
     *
     * The Node class bundles a number of functions for bus connections,
     * printing, names, etc.
     */
    template <typename T>
    class Node
    {
    public:


        //! \brief Constructs a new, unitialized node
        Node(): name_(""), busIndex_(-1) {}


        virtual ~Node() {}


        //! \brief Connects the load to a bus
        T& bus(Bus const& bus)
        {
            busIndex_ = bus.index;
            return dynamic_cast<T&>(*this);
        }


        std::size_t busIndex() const
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
        T& name(std::string const& name)
        {
            name_ = name;
            return dynamic_cast<T&>(*this);
        }


    private:


        //! \brief The name, for readability
        std::string name_;


        //! \brief Index of the bus we're connected to.
        std::ptrdiff_t busIndex_;
    };
} // namespace fPotencia

#endif // FPOTENCIA_NODE_H
