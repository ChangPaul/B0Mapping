/*
 * sphaFns.h
 * Definitions of the spherical harmonic functions
 *
 *  Created on: Jan 20, 2016
 *      Author: Paul Chang
 */

#ifndef SPHA_H_INCLUDED
#define SPHA_H_INCLUDED
 
#define NUM_SPHA 64

typedef double ( *basisFn )( double x, double y, double z );

inline double rho2               (double x, double y)         { return (x*x + y*y); }
inline double rho4               (double x, double y)         { return rho2(x,y)*rho2(x,y); }
inline double rho6               (double x, double y)         { return rho4(x,y)*rho2(x,y); }
inline double rho8               (double x, double y)         { return rho4(x,y)*rho4(x,y); }

inline double ordZero            (double x, double y, double z){ return 1; }
inline double ordOne_degNOne     (double x, double y, double z){ return y; }
inline double ordOne_degZero     (double x, double y, double z){ return z; }
inline double ordOne_degOne      (double x, double y, double z){ return x; }

inline double ordTwo_degNTwo     (double x, double y, double z){ return x*y; }
inline double ordTwo_degNOne     (double x, double y, double z){ return y*z; }
inline double ordTwo_degZero     (double x, double y, double z){ return (z*z - 0.5*rho2(x,y)); }
inline double ordTwo_degOne      (double x, double y, double z){ return x*z; }
inline double ordTwo_degTwo      (double x, double y, double z){ return (x*x - y*y); }

inline double ordThree_degNThree (double x, double y, double z){ return (3*x*x - y*y)*y; }
inline double ordThree_degNTwo   (double x, double y, double z){ return x*y*z; }
inline double ordThree_degNOne   (double x, double y, double z){ return y*(z*z - 0.25*rho2(x,y)); }
inline double ordThree_degZero   (double x, double y, double z){ return z*(z*z - 1.50*rho2(x,y)); }
inline double ordThree_degOne    (double x, double y, double z){ return x*(z*z - 0.25*rho2(x,y)); }
inline double ordThree_degTwo    (double x, double y, double z){ return (x*x - y*y)*z; }
inline double ordThree_degThree  (double x, double y, double z){ return (x*x - 3*y*y)*x; }

inline double ordFour_degNFour   (double x, double y, double z){ return (  x*x - y*y)*x*y; }
inline double ordFour_degNThree  (double x, double y, double z){ return (3*x*x - y*y)*y*z; }
inline double ordFour_degNTwo    (double x, double y, double z){ return x*y*        (z*z - 1.0/6*rho2(x,y)); }
inline double ordFour_degNOne    (double x, double y, double z){ return y*z*        (z*z - 3.0/4*rho2(x,y)); }
inline double ordFour_degZero    (double x, double y, double z){ return z*z*        (z*z -   3*rho2(x,y)) + 3.0/8*rho4(x,y); }
inline double ordFour_degOne     (double x, double y, double z){ return x*z*        (z*z - 3.0/4*rho2(x,y)); }
inline double ordFour_degTwo     (double x, double y, double z){ return (x*x - y*y)*(z*z - 1.0/6*rho2(x,y)); }
inline double ordFour_degThree   (double x, double y, double z){ return (x*x - 3*y*y)*x*z; }
inline double ordFour_degFour    (double x, double y, double z){ return x*x*(x*x - 6*y*y) + y*y*y*y; }

inline double ordFive_degNFive   (double x, double y, double z){ return y*(y*y*(y*y - 10*x*x) + 5*x*x*x*x); }
inline double ordFive_degNFour   (double x, double y, double z){ return z*(  x*x - y*y)*x*y; }
inline double ordFive_degNThree  (double x, double y, double z){ return y*(3*x*x - y*y)*(z*z - 1.0/8*rho2(x,y)); }
inline double ordFive_degNTwo    (double x, double y, double z){ return z*x*y*        (z*z - 0.5*rho2(x,y)); }
inline double ordFive_degNOne    (double x, double y, double z){ return y* (z*z* (z*z - 1.5*rho2(x,y)) +  1.0/8*rho4(x,y)); }
inline double ordFive_degZero    (double x, double y, double z){ return z* (z*z* (z*z - 5.0*rho2(x,y)) + 15.0/8*rho4(x,y)); }
inline double ordFive_degOne     (double x, double y, double z){ return x* (z*z* (z*z - 1.5*rho2(x,y)) +  1.0/8*rho4(x,y)); }
inline double ordFive_degTwo     (double x, double y, double z){ return z*(x*x - y*y)*(z*z - 0.5*rho2(x,y)); }
inline double ordFive_degThree   (double x, double y, double z){ return x*(x*x - 3*y*y)*(z*z - 1.0/8*rho2(x,y)); }
inline double ordFive_degFour    (double x, double y, double z){ return z*(x*x*(x*x - 6*y*y) + y*y*y*y); }
inline double ordFive_degFive    (double x, double y, double z){ return x*(x*x*(x*x - 10*y*y) + 5*y*y*y*y); }

inline double ordSix_degNSix    (double x, double y, double z){ return x*y*(x*x - 3.0*y*y)*(3.0*x*x - y*y); }
inline double ordSix_degNFive   (double x, double y, double z){ return ordFive_degNFive(x,y,z)*z; }
inline double ordSix_degNFour   (double x, double y, double z){ return ordFour_degNFour(x,y,z)*(z*z - 1.0/10*rho2(x,y)); }
inline double ordSix_degNThree  (double x, double y, double z){ return y*z*(3*x*x - y*y)*(z*z - 3.0/8*rho2(x,y)); }
inline double ordSix_degNTwo    (double x, double y, double z){ return x*y*(z*z* (z*z - rho2(x,y)) + 1.0/16*rho4(x,y)); }
inline double ordSix_degNOne    (double x, double y, double z){ return y*z* (z*z* (z*z -  5.0/2*rho2(x,y)) +  5.0/8*rho4(x,y)); }
inline double ordSix_degZero    (double x, double y, double z){ return z*z* (z*z* (z*z - 15.0/2*rho2(x,y)) + 45.0/8*rho4(x,y)) - 5.0/16*rho2(x,y)*rho4(x,y); }
inline double ordSix_degOne     (double x, double y, double z){ return x*z* (z*z* (z*z -  5.0/2*rho2(x,y)) +  5.0/8*rho4(x,y)); }
inline double ordSix_degTwo     (double x, double y, double z){ return (x*x - y*y)*(z*z* (z*z - rho2(x,y)) + 1.0/16*rho4(x,y)); }
inline double ordSix_degThree   (double x, double y, double z){ return x*z*(x*x - 3.0*y*y)*(z*z - 3.0/8*rho2(x,y)); }
inline double ordSix_degFour    (double x, double y, double z){ return ordFour_degFour(x,y,z)*(z*z - 1.0/10*rho2(x,y)); }
inline double ordSix_degFive    (double x, double y, double z){ return ordFive_degFive(x,y,z)*z; }
inline double ordSix_degSix     (double x, double y, double z){ return (x*x - y*y)*(x*x*(x*x - 14.0*y*y) + y*y*y*y); }

inline double ordSeven_degNSeven(double x, double y, double z){ return (y*(y*y*y*y*(y*y - 21*x*x) + (35*y*y - 7*x*x)*x*x*x*x)); }
inline double ordSeven_degNSix  (double x, double y, double z){ return ordSix_degNSix(x,y,z)*z; }
inline double ordSeven_degNFive (double x, double y, double z){ return ordFive_degNFive(x,y,z)*(z*z - 1.0/12*rho2(x,y)); }
inline double ordSeven_degNFour (double x, double y, double z){ return ordFive_degNFour(x,y,z)*(z*z - 3.0/10*rho2(x,y)); }
inline double ordSeven_degNThree(double x, double y, double z){ return ordThree_degNThree(x,y,z)*(z*z*(z*z - 3.0/4*rho2(x,y)) + 3.0/80*rho4(x,y)); }
inline double ordSeven_degNTwo  (double x, double y, double z){ return z*x*y*(z*z*(z*z - 5.0/3*rho2(x,y)) + 5.0/16*rho4(x,y)); }
inline double ordSeven_degNOne  (double x, double y, double z){ return y*(z*z*z*z*(z*z - 15.0/4*rho2(x,y)) + 15.0/8*z*z*rho4(x,y) - 5.0/64*rho6(x,y)); }
inline double ordSeven_degZero  (double x, double y, double z){ return z*(z*z*z*z*(z*z - 21.0/2*rho2(x,y)) + 105.0/8*z*z*rho4(x,y) - 35.0/16*rho6(x,y)); }
inline double ordSeven_degOne   (double x, double y, double z){ return x*(z*z*z*z*(z*z - 15.0/4*rho2(x,y)) + 15.0/8*z*z*rho4(x,y) - 5.0/64*rho6(x,y)); }
inline double ordSeven_degTwo   (double x, double y, double z){ return z*(x*x - y*y)*(z*z*(z*z - 5.0/3*rho2(x,y)) + 5/16*rho4(x,y)); }
inline double ordSeven_degThree (double x, double y, double z){ return ordThree_degThree(x,y,z)*(z*z*(z*z - 3.0/4*rho2(x,y)) + 3.0/80*rho4(x,y)); }
inline double ordSeven_degFour  (double x, double y, double z){ return ordFive_degFour(x,y,z)*(z*z - 3.0/10*rho2(x,y)); }
inline double ordSeven_degFive  (double x, double y, double z){ return ordFive_degFive(x,y,z)*(z*z - 1.0/12*rho2(x,y)); }
inline double ordSeven_degSix   (double x, double y, double z){ return ordSix_degSix(x,y,z)*z; }
inline double ordSeven_degSeven (double x, double y, double z){ return x*(x*x*x*x*(x*x - 21*y*y) + (35*x*x - 7*y*y)*y*y*y*y); }

const basisFn sphFn [NUM_SPHA] = {
	ordZero, ordOne_degZero, ordOne_degOne, ordOne_degNOne,
    ordTwo_degZero, ordTwo_degOne, ordTwo_degNOne, ordTwo_degTwo, ordTwo_degNTwo,
	ordThree_degZero, ordThree_degOne, ordThree_degNOne, ordThree_degTwo, ordThree_degNTwo, ordThree_degThree, ordThree_degNThree,
	ordFour_degZero, ordFour_degOne, ordFour_degNOne, ordFour_degTwo, ordFour_degNTwo,
	ordFour_degThree, ordFour_degNThree, ordFour_degFour, ordFour_degNFour,
	ordFive_degZero, ordFive_degOne, ordFive_degNOne, ordFive_degTwo, ordFive_degNTwo,
	ordFive_degThree, ordFive_degNThree, ordFive_degFour, ordFive_degNFour,
	ordFive_degFive, ordFive_degNFive,
	ordSix_degZero, ordSix_degOne, ordSix_degNOne, ordSix_degTwo, ordSix_degNTwo,
	ordSix_degThree, ordSix_degNThree, ordSix_degFour, ordSix_degNFour,
	ordSix_degFive, ordSix_degNFive, ordSix_degSix, ordSix_degNSix,
	ordSeven_degZero, ordSeven_degOne, ordSeven_degNOne, ordSeven_degTwo, ordSeven_degNTwo,
	ordSeven_degThree, ordSeven_degNThree, ordSeven_degFour, ordSeven_degNFour,
	ordSeven_degFive, ordSeven_degNFive, ordSeven_degSix, ordSeven_degNSix,
	ordSeven_degSeven, ordSeven_degNSeven
};

#endif