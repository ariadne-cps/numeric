/***************************************************************************
 *            numeric/api.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/api.hpp
 *  @brief Export and visibility macros for the Numeric library.
 */

#ifndef ARIADNE_NUMERIC_API_HPP
#define ARIADNE_NUMERIC_API_HPP

#if defined(_WIN32)
#if defined(ARIADNE_NUMERIC_STATIC) || defined(ARIADNE_NUMERIC_EMBEDDED)
#define ARIADNE_NUMERIC_API
#elif defined(ARIADNE_NUMERIC_BUILD)
#define ARIADNE_NUMERIC_API __declspec(dllexport)
#else
#define ARIADNE_NUMERIC_API __declspec(dllimport)
#endif
#else
#define ARIADNE_NUMERIC_API
#endif

#endif
