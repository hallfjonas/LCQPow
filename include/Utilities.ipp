/*
 *	This file is an adaptation of the implementation in qpOASES (Utils.ipp).
 */

/*
 *	This file is part of lcqpOASES.
 *
 *	lcqpOASES -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	lcqpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	lcqpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with lcqpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Types.hpp"
#include "Utilities.hpp"

#include <math.h>

namespace lcqpOASES {
    /*
    *   i s E q u a l
    */
    inline bool isEqual(real_t x, real_t y, real_t TOL)
    {
        if ( getAbs(x-y) <= TOL )
            return true;
        else
            return false;
    }


    /*
    *   i s E q u a l
    */
    inline bool isEqual(real_t x, real_t y)
    {
        return isEqual(x, y, Utilities::ZERO);
    }


    /*
    *   i s Z e r o
    */
    inline bool isZero(real_t x, real_t TOL)
    {
        if ( getAbs(x) <= TOL )
            return true;
        else
            return false;
    }


    /*
    *   i s Z e r o
    */
    inline bool isZero(real_t x)
    {
        return isZero(x, Utilities::ZERO);
    }

    /*
    *   g e t S i g n
    */
    inline real_t getSign(real_t arg)
    {
        if ( arg >= 0.0 )
            return 1.0;
        else
            return -1.0;
    }


    /*
    *   g e t M a x
    */
    inline int_t getMax(int_t x, int_t y)
    {
        return (y<x) ? x : y;
    }


    /*
    *   g e t M i n
    */
    inline int_t getMin(int_t x, int_t y)
    {
        return (y>x) ? x : y;
    }



    /*
    *   g e t M a x
    */
    inline real_t getMax(real_t x, real_t y)
    {
        return (y<x) ? x : y;
    }


    /*
    *   g e t M i n
    */
    inline real_t getMin(real_t x, real_t y)
    {
        return (y>x) ? x : y;
    }


    /*
    *   g e t A b s
    */
    inline real_t getAbs(real_t x)
    {
        #ifdef __NO_FMATH__
        return (x>=0.0) ? x : -x;
        #else
        return fabs(x);
        #endif
    }


    /*
    *   g e t S q r t
    */
    inline real_t getSqrt(real_t x)
    {
        #ifdef __NO_FMATH__
        return sqrt(x); /* put your custom sqrt-replacement here */
        #else
        return sqrt(x);
        #endif
    }
}


/*
 *	end of file
 */
