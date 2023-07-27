/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>

#include "FiniteDiffDeriv.H"

// five point stencil first and second derivatives using arbitrary spacing
//
// dxN is the distance between the points with function values yN and y0
//
double FivePointStencilFirstDeriv( const double dx1,
                                   const double dx2,
                                   const double dx3,
                                   const double dx4,
                                   const double y0,
                                   const double y1,
                                   const double y2,
                                   const double y3,
                                   const double y4 )
{
  return
    -(dx1*dx2*dx3+dx1*dx2*dx4+dx1*dx3*dx4+dx2*dx3*dx4)/dx1/dx2/dx3/dx4*y0
    -dx2*dx3*dx4/dx1/(dx1-dx2)/(dx1-dx3)/(dx1-dx4)*y1
    +dx1*dx3*dx4/(dx1-dx2)/dx2/(dx2-dx3)/(dx2-dx4)*y2
    -dx1*dx2*dx4/dx3/(dx3-dx1)/(dx3-dx2)/(dx3-dx4)*y3
    +dx1*dx2*dx3/(dx1-dx4)/(dx2-dx4)/(dx3-dx4)/dx4*y4;
}

double FivePointStencilSecondDeriv( const double dx1,
                                    const double dx2,
                                    const double dx3,
                                    const double dx4,
                                    const double y0,
                                    const double y1,
                                    const double y2,
                                    const double y3,
                                    const double y4 )
{
  return
    2.*((dx1*dx2+dx1*dx3+dx2*dx3+dx1*dx4+dx2*dx4+dx3*dx4)/dx1/dx2/dx3/dx4*y0
        +(dx2*dx3+dx2*dx4+dx3*dx4)/dx1/(dx1-dx2)/(dx1-dx3)/(dx1-dx4)*y1
        -(dx1*dx3+dx1*dx4+dx3*dx4)/(dx1-dx2)/dx2/(dx2-dx3)/(dx2-dx4)*y2
        +(dx1*dx2+dx1*dx4+dx2*dx4)/dx3/(dx3-dx1)/(dx3-dx2)/(dx3-dx4)*y3
        -(dx1*dx2+dx1*dx3+dx2*dx3)/(dx1-dx4)/(dx2-dx4)/(dx3-dx4)/dx4*y4);
}

// seven point stencil derivatives using equal spacing:
//
// a0 f(x+3a)+a1 f(x+2a)+a2 f(x+a) +a3 f(x)   +a4 f(x-a) +a5 f(x-2a)+a6 f(x-3a)
//
//                 a0 a1  a2 a3 a4 a5 a6
//  60 a   f'   =>  1  -9  45    0 -45   9 -1
// 180 a^2 f''  =>  2 -27 270 -490 270 -27  2
//   8 a^3 f''' => -1   8 -13    0  13  -8  1
//
// a0 f(x+2a)+a1 f(x+a) +a2 f(x)   +a3 f(x-a) +a4 f(x-2a)+a5 f(x-3a)+a6 f(x-4a)
//
//                 a0 a1  a2 a3 a4 a5 a6
//  60 a   f'   =>  -2  24   35 -80  30  -8  1
// 180 a^2 f''  => -13 228 -420 200  15 -12  2
//   8 a^3 f''' =>   1   8  -35  48 -29   8 -1
//
// a0 f(x+a) +a1 f(x)   +a2 f(x-a) +a3 f(x-2a)+a4 f(x-3a)+a5 f(x-4a)+a6 f(x-5a)
//
//                 a0 a1  a2 a3 a4 a5 a6
//  60 a   f'   =>  10   77 -150 100  -50 15  -2
// 180 a^2 f''  => 137 -147 -255 470 -285 93 -13
//   8 a^3 f''' =>  15  -56   83 -64   29 -8   1
//
// a0 f(x)   +a1 f(x-a) +a2 f(x-2a)+a3 f(x-3a)+a4 f(x-4a)+a5 f(x-5a)+a6 f(x-6a)
//
//                 a0 a1  a2 a3 a4 a5 a6
//  60 a   f'   => 147  -360  450  -400  225  -72  10
// 180 a^2 f''  => 812 -3132 5265 -5080 2970 -972 137
//   8 a^3 f''' =>  49  -232  461  -496  307 -104  15
//
// f is a vector of length 7 containing function values at an increment of a
// zero is the index corresponding to f(x)
// order is the derivative order
//
double SevenPointStencilDeriv( std::vector<double> & f,
                               double a,
                               int zero,
                               int order )
{
  if (zero < 0 || zero > 3)
    throw std::logic_error("bad zero index in SevenPointStencilDeriv");

  if (order < 1 || order > 3)
    throw std::logic_error("bad order in SevenPointStencilDeriv");

  if (f.size() < 7)
    throw std::logic_error("function vector too small in SevenPointStencilDeriv");

  if (zero == 3) {
    if (order == 1)
      return (f[0]-9.*f[1]+45.*f[2]-45.*f[4]+9.*f[5]-f[6])/a/60.;
    else if (order == 2)
      return (2.*f[0]-27.*f[1]+270.*f[2]-490.*f[3]+270.*f[4]-27.*f[5]+2.*f[6])/a/a/180.;
    else
      return (-f[0]+8.*f[1]-13.*f[2]+13.*f[4]-8.*f[5]+f[6])/a/a/a/8.;
  }
  else if (zero == 2) {
    if (order == 1)
      return (-2.*f[0]+24.*f[1]+35.*f[2]-80.*f[3]+30.*f[4]-8.*f[5]+f[6])/a/60.;
    else if (order == 2)
      return (-13.*f[0]+228.*f[1]-420.*f[2]+200.*f[3]+15.*f[4]-12.*f[5]+2.*f[6])/a/a/180.;
    else
      return (f[0]+8.*f[1]-35.*f[2]+48.*f[3]-29.*f[4]+8.*f[5]-f[6])/a/a/a/8.;
  }
  else if (zero == 1) {
    if (order == 1)
      return (10.*f[0]+77.*f[1]-150.*f[2]+100.*f[3]-50.*f[4]+15.*f[5]-2.*f[6])/a/60.;
    else if (order == 2)
      return (137.*f[0]-147.*f[1]-255.*f[2]+470.*f[3]-285.*f[4]+93.*f[5]-13.*f[6])/a/a/180.;
    else
      return (15.*f[0]-56.*f[1]+83.*f[2]-64.*f[3]+29.*f[4]-8.*f[5]+f[6])/a/a/a/8.;
  }
  else {
    if (order == 1)
      return (147.*f[0]-360.*f[1]+450.*f[2]-400.*f[3]+225.*f[4]-72.*f[5]+10.*f[6])/a/60.;
    else if (order == 2)
      return (812.*f[0]-3132.*f[1]+5265.*f[2]-5080.*f[3]+2970.*f[4]-972.*f[5]+137.*f[6])/a/a/180.;
    else
      return (49.*f[0]-232.*f[1]+461.*f[2]-496.*f[3]+307.*f[4]-104.*f[5]+15.*f[6])/a/a/a/8.;
  }

}

// five point stencil derivatives using equal spacing:
//
// a0 f(x+2a) + a1 f(x+a) + a2 f(x)    + a3 f(x-a)  + a4 f(x-2a)
//                a0 a1  a2 a3 a4
// 12 a   f'   => -1  8   0 -8  1
// 12 a^2 f''  => -1 16 -30 16 -1
//  2 a^3 f''' =>  1 -2   0  2 -1
//
// a0 f(x+a)  + a1 f(x)   + a2 f(x-a)  + a3 f(x-2a) + a4 f(x-3a)
//                a0  a1  a2 a3 a4
// 12 a   f'   =>  3  10 -18  6 -1
// 12 a^2 f''  => 11 -20   6  4 -1
//  2 a^3 f''' =>  3 -10  12 -6  1
//
// a0 f(x)    + a1 f(x-a) + a2 f(x-2a) + a3 f(x-3a) + a4 f(x-4a)
//                a0   a1  a2  a3 a4
// 12 a   f'   => 25  -48  36 -16  3
// 12 a^2 f''  => 35 -104 114 -56 11
//  2 a^3 f''' =>  5  -18  24 -14  3
//
// f is a vector of length 5 containing function values at an increment of a
// zero is the index corresponding to f(x)
// order is the derivative order
//
double FivePointStencilDeriv( std::vector<double> & f,
                              double a,
                              int zero,
                              int order )
{
  if (zero < 0 || zero > 2)
    throw std::logic_error("bad zero index in FivePointStencilDeriv");

  if (order < 1 || order > 3)
    throw std::logic_error("bad order in FivePointStencilDeriv");

  if (f.size() < 5)
    throw std::logic_error("function vector too small in FivePointStencilDeriv");

  if (zero == 2) {
    if (order == 1)
      return (-f[0]+8.*f[1]-8.*f[3]+f[4])/a/12.;
    else if (order == 2)
      return (-f[0]+16.*f[1]-30.*f[2]+16.*f[3]-f[4])/a/a/12.;
    else
      return ( f[0]-2.*f[1]+2.*f[3]-f[4])/a/a/a/2.;
  }
  else if (zero == 1) {
    if (order == 1)
      return (3.*f[0]+10.*f[1]-18.*f[2]+6.*f[3]-f[4])/a/12.;
    else if (order == 2)
      return (11.*f[0]-20.*f[1]+6.*f[2]+4.*f[3]-f[4])/a/a/12.;
    else
      return (3.*f[0]-10.*f[1]+12.*f[2]-6.*f[3]+f[4])/a/a/a/2.;
  }
  else {
    if (order == 1)
      return (25.*f[0]-48.*f[1]+36.*f[2]-16.*f[3]+3.*f[4])/a/12.;
    else if (order == 2)
      return (35.*f[0]-104.*f[1]+114.*f[2]-56.*f[3]+11.*f[4])/a/a/12.;
    else
      return (5.*f[0]-18.*f[1]+24.*f[2]-14.*f[3]+3.*f[4])/a/a/a/2.;
  }

}
