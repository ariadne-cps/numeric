/***************************************************************************
 *            numeric_module.cpp
 *
 *  Copyright  2007-26  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */

#include "pybind11.hpp"

void numeric_submodule(pybind11::module& module);

void numeric_module(pybind11::module& module) {
    numeric_submodule(module);
}
