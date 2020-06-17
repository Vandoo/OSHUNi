///////////////////////////////////////////////////////////
//   
//   This file contains the definitions for the  
//   classes ElectricField, MagneticField, Current, 
//   Spatial_Advection and MaxwellEq.
///////////////////////////////////////////////////////////
//
//   classes :
//
//   1. class ElectricField :
//      The class ElectricField  stores a reference to the slope Yh 
//      and steps it calls Exyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      electric field on the distribution function.
//
//   2. class MagneticField :
//      The class Magneticfield stores a reference to the slope Yh 
//      and steps it calls Bxyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      magnetic field on the distribution function.
//
//   3. class Current :
//      The class Current stores a reference to the slope Yh and
//      in subsequent steps it calls Jx(Yin), where Yin is the 
//      state, in order  to calculate effect of the current
//      on the Electric field.
//
//   4. Spatial_Advection :
//      The class Spatial_Advection  stores a reference to the slope Yh and
//      in subsequent steps it calls Axy(Yin), where Yin is the 
//      state, in order to calculate effect of the spatial 
//      advection on the distribution function. 
//
//   5. class MaxwellEq :
//      The class MaxwellEq stores a reference to the slope Yh and
//      in subsequent steps it calls Ampere(Yin) and Faraday(Yin),
//      where Yin is the state, in order to update the electromagnetic
//      fields 
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard Libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <math.h>
    #include <float.h>

//  My Libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-vlasovmax.h"




//**************************************************************
//**************************************************************
//   Definition for the Electric Field in the x direction
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    ElectricField::ElectricField(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(l0+1,m0+1), A2(l0+1,m0+1),
         B1(l0+1),      B2(l0+1),
         C1(l0+1),      C3(l0+1),
         C2(l0+1,m0+1), C4(l0+1,m0+1),
         H(Yslope.SH(0,0)), G(Yslope.SH(0,0)),TMP(Yslope.SH(0,0)),
		 TMP2(Yslope.SH(0,0)),
         Hp0(l0+1),
         pr(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         invpr(pr)
         {
         complex<double> lc, mc, tmp;

//       Calculate the "A1, A2" parameters                      //WHY- Eq.(38)
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
			 tmp = 1.0 /(2.0*lc + 1.0);
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 mc = static_cast< complex<double> >(m);
                 A1(l,m) = -(lc + 1.0 - mc) * tmp;
                 A2(l,m) = -(lc + mc)       * tmp;
				 // Calculate the "C2, C4" parameters
				 if(l>1 && m>1){
					C2(l,m) = (lc-mc+2.0)*(lc-mc+1.0) * tmp;
					C4(l,m) =-(lc+mc-1.0)*(lc+mc)     * tmp;
				 }
             }
			 // Calculate the "C1, C3" parameters
             C3[l] = (0.25/pr.dx()) * tmp;
			 C1[l] = -C3[l];
			 // Calculate the "B1, B2" parameters
             B1[l] =  C3[l]*2.0 *lc * (lc+1.0);
             B2[l] = -B1[l];
		 }
		 A1(0,0) = -1.0;
         A1 *= 0.5 / pr.dx(); A2 *= 0.5 / pr.dx();
		 C2 *= 0.25/ pr.dx(); C4 *= 0.25/ pr.dx();
         C1[0] = -0.25/pr.dx();
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "H" parameters at the p0 cell            //WHY- (34)
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Hp0 = 1.0 / pr(0); 
         for (size_t l(1); l < l0+1; ++l) { Hp0[l] = Hp0[l-1]*(pr(0)/pr(1)) * (2.0*l+1.0)/(2.0*l-1);}
         Hp0 *= -2.0 * pr.dx(); 
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Create inverted momentum axis   
         for (size_t i(0); i<invpr.dim(); ++i) { invpr(i) = 1.0/pr(i); } 

		 for (size_t l(0); l<=l0; ++l){
			valarray< complex<double> > pl(Inputdata::IN().inp().pr.dim());
			valarray< complex<double> > inv_pl(Inputdata::IN().inp().pr.dim());
			for (size_t i(0); i<pr.dim(); ++i) {
			  if(l == 0)
				  pl[i] = pr(i);
			  else
				  pl[i] = pow(pr(i),l+1);
			  inv_pl[i] = 1.0 / pl[i];
			}
			pls.push_back(pl);
			inv_pls.push_back(inv_pl);
		 }
		 // output for tests
		 /*for(size_t i(0); i<5;++i)
			for(size_t j(0); j<4;++j){
				cout<<i<<": "<<pls[i][j]<<" "<<inv_pls[i][j]<<endl;
		 }*/
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Ex(Stat& Yin, Matrix2D< complex<double> >& Ex){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0);             Yh.SH(1,0) += G.mxy_matrix(Ex);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);   Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0);   Yh.SH(2,0) += G.mxy_matrix(Ex);

//      m = 0, l = 2
        MakeGH(2,Yin.SH(2,0));
        Ex *= A2(2,0) / A1(1,0);   Yh.SH(1,0) += H.mxy_matrix(Ex);

//      m = 0, l = 3
        MakeGH(3,Yin.SH(3,0));
        Ex *= A2(3,0) / A2(2,0);   Yh.SH(2,0) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ex *= A1(1,1) / A2(3,0);   Yh.SH(2,1) += G.mxy_matrix(Ex); 

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ex *= A2(2,1) / A1(1,1);   Yh.SH(1,1) += H.mxy_matrix(Ex);

//      m = 1, l = 3
        MakeGH(3,Yin.SH(3,1));
        Ex *= A2(3,1) / A2(2,1);   Yh.SH(2,1) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Ex1D(Stat& Yin, Matrix2D< complex<double> >& Ex){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0);             Yh.SH(1,0) += G.mxy_matrix(Ex);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);   Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0);   Yh.SH(2,0) += G.mxy_matrix(Ex);

//      m = 0, l = 2
        MakeGH(2,Yin.SH(2,0));
        Ex *= A2(2,0) / A1(1,0);   Yh.SH(1,0) += H.mxy_matrix(Ex);

//      m = 0, l = 3
        MakeGH(3,Yin.SH(3,0));
        Ex *= A2(3,0) / A2(2,0);   Yh.SH(2,0) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Eyz(Stat& Yin, Matrix2D< complex<double> >& Em, 
                                                 Matrix2D< complex<double> >& Ep){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Em *= C1[0];            Yh.SH(1,1) += G.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Em *= C1[1]   / C1[0];  Yh.SH(2,1) += G.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        MakeGH(2,Yin.SH(2,0));
        Em *= C3[2]   / C1[1];  Yh.SH(1,1) += H.mxy_matrix(Em);

//      m = 0,  l = 3 
        MakeGH(3,Yin.SH(3,0));
        Em *= C3[3]   / C3[2];  Yh.SH(2,1) += H.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ep *= B2[1];             H = H.mxy_matrix(Ep); Yh.SH(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;              Yh.SH(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];     G = G.mxy_matrix(Ep); Yh.SH(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Yh.SH(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(3,Yin.SH(3,1));
        Em *= C3[3]   / C1[1];   TMP = H;              Yh.SH(2,2) += TMP.mxy_matrix(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxy_matrix(Ep); Yh.SH(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(2,Yin.SH(2,2));
        Ep *= C4(2,2) / B2[3];        Yh.SH(1,1) += H.mxy_matrix(Ep); 

//      m = 2, l = 3
        MakeGH(3,Yin.SH(3,2));
        Ep *= C4(3,2) / C4(2,2);      Yh.SH(2,1)  += H.mxy_matrix(Ep);

        if ( m0 > 2) { 
//          m = 3, l = 3
            MakeGH(3,Yin.SH(3,3));
            Ep *= C4(3,3) / C4(3,2);  Yh.SH(2,2)  += H.mxy_matrix(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Exyz1D(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        //Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix());
        //Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *= (-1.0)*ii; Em += Yin.EMF().Ey().matrix();
        //Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *=  ii;       Ep += Yin.EMF().Ey().matrix();
		// --- e-i
        Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix()); Ex *= (-1.0);
        Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *= ii;        Em -= Yin.EMF().Ey().matrix();
        Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *= (-1.0)*ii; Ep -= Yin.EMF().Ey().matrix();


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0); TMP = G; Yh.SH(1,0) += G.mxy_matrix(Ex);
//        Em *= C1[0];            Yh.SH(1,1) += TMP.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);          Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0); TMP = G; Yh.SH(2,0) += G.mxy_matrix(Ex);
//        Em *= C1[1]   / C1[0]  ;          Yh.SH(2,1) += TMP.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        for (size_t l(2); l < l0; ++l){
            MakeGH(l,Yin.SH(l,0));
            Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Yh.SH(l-1,0) += H.mxy_matrix(Ex);
//            Em *= C3[l]   / C1[l-1];             Yh.SH(l-1,1) += TMP.mxy_matrix(Em);
            Ex *= A1(l,0) / A2(l,0);   TMP = G;  Yh.SH(l+1,0) += G.mxy_matrix(Ex);
//            Em *= C1[l]   / C3[l];               Yh.SH(l+1,1) += TMP.mxy_matrix(Em);
        }
//      m = 0,  l = l0
        MakeGH(l0,Yin.SH(l0,0));
        Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Yh.SH(l0-1,0) += H.mxy_matrix(Ex);
//        Em *= C3[l0]   / C1[l0-1];             Yh.SH(l0-1,1) += TMP.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// This is the end for the 1D version
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& ElectricField::Exyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        //Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix());
        //Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *= (-1.0)*ii; Em += Yin.EMF().Ey().matrix();
        //Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *=  ii;       Ep += Yin.EMF().Ey().matrix();
		// --- e-i
        Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix()); Ex *= (-1.0);
        Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *=        ii; Em -= Yin.EMF().Ey().matrix();
        Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *= (-1.0)*ii; Ep -= Yin.EMF().Ey().matrix();


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0); TMP = G; Yh.SH(1,0) += G.mxy_matrix(Ex);
        Em *= C1[0];            Yh.SH(1,1) += TMP.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);          Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0); TMP = G; Yh.SH(2,0) += G.mxy_matrix(Ex);
        Em *= C1[1]   / C1[0]  ;          Yh.SH(2,1) += TMP.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        for (size_t l(2); l < l0; ++l){
            MakeGH(l,Yin.SH(l,0));
            Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Yh.SH(l-1,0) += H.mxy_matrix(Ex);
            Em *= C3[l]   / C1[l-1];             Yh.SH(l-1,1) += TMP.mxy_matrix(Em);
            Ex *= A1(l,0) / A2(l,0);   TMP = G;  Yh.SH(l+1,0) += G.mxy_matrix(Ex);
            Em *= C1[l]   / C3[l];               Yh.SH(l+1,1) += TMP.mxy_matrix(Em);
        }
//      m = 0,  l = l0
        MakeGH(l0,Yin.SH(l0,0));
        Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Yh.SH(l0-1,0) += H.mxy_matrix(Ex);
        Em *= C3[l0]   / C1[l0-1];             Yh.SH(l0-1,1) += TMP.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ep *= B2[1];              H = H.mxy_matrix(Ep); Yh.SH(0,0) += H.Re();
        Ex *= A1(1,1) / A2(l0,0); TMP = G;              Yh.SH(2,1) += G.mxy_matrix(Ex); 
        Em *= C1[1] / C3[l0];     G = TMP;              Yh.SH(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];      G = G.mxy_matrix(Ep); Yh.SH(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ex *= A2(2,1) / A1(1,1); TMP = H;              Yh.SH(1,1) += TMP.mxy_matrix(Ex);
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Yh.SH(1,0) += H.Re();
        Ex *= A1(2,1) / A2(2,1); TMP = G;              Yh.SH(3,1) += G.mxy_matrix(Ex);
        Em *= C1[2]   / C1[1];   G = TMP;              Yh.SH(3,2) += TMP.mxy_matrix(Em);
        Ep *= B1[2]   / B2[2];   G = G.mxy_matrix(Ep); Yh.SH(3,0) += G.Re();

//      m = 1, 1 < l < l0
        for (size_t l(3); l < l0; ++l){
            MakeGH(l,Yin.SH(l,1));
            Ex *= A2(l,1) / A1(l-1,1); TMP = H;              Yh.SH(l-1,1) += H.mxy_matrix(Ex);
            Em *= C3[l]   / C1[l-1];   H = TMP;              Yh.SH(l-1,2) += TMP.mxy_matrix(Em);
            Ep *= B2[l]   / B1[l-1];   H = H.mxy_matrix(Ep); Yh.SH(l-1,0) += H.Re();
            Ex *= A1(l,1) / A2(l,1);   TMP = G;              Yh.SH(l+1,1) += G.mxy_matrix(Ex);
            Em *= C1[l]   / C3[l];     G = TMP;              Yh.SH(l+1,2) += TMP.mxy_matrix(Em);
            Ep *= B1[l]   / B2[l];     G = G.mxy_matrix(Ep); Yh.SH(l+1,0) += G.Re();
         }
//       m = 1,  l = l0
         MakeGH(l0,Yin.SH(l0,1));
         Ex *= A2(l0,1) / A1(l0-1,1); TMP = H;              Yh.SH(l0-1,1) += H.mxy_matrix(Ex);
         Em *= C3[l0]   / C1[l0-1];   H = TMP;              Yh.SH(l0-1,2) += TMP.mxy_matrix(Em);
         Ep *= B2[l0]   / B1[l0-1];   H = H.mxy_matrix(Ep); Yh.SH(l0-1,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        C4(l0,1) = B2[l0];
        for (size_t m(2); m < m0; ++m){
//          m > 1 , l = m
            MakeGH(m,Yin.SH(m,m));
            Ep *= C4(m,m) / C4(l0,m-1);              Yh.SH(m-1,m-1) += H.mxy_matrix(Ep); 
            Ex *= A1(m,m) / A2(l0,m-1); TMP = G;     Yh.SH(m+1,m  ) += G.mxy_matrix(Ex); 
            Em *= C1[m]   / C3[l0];     G = TMP;     Yh.SH(m+1,m+1) += TMP.mxy_matrix(Em); 
            Ep *= C2(m,m) / C4(m,m);                 Yh.SH(m+1,m-1) += G.mxy_matrix(Ep);

//          m > 1 , l = m+1
            MakeGH(m+1,Yin.SH(m+1,m));
            Ex *= A2(m+1,m) / A1(m,m);   TMP = H;     Yh.SH(m  ,m  ) += TMP.mxy_matrix(Ex);
            Ep *= C4(m+1,m) / C2(m,m);                Yh.SH(m  ,m-1) += H.mxy_matrix(Ep);
            if ( m+1 < l0) { //always true except when m = m0-1 = l0-1 
                Ex *= A1(m+1,m) / A2(m+1,m); TMP = G;     Yh.SH(m+2,m  ) += G.mxy_matrix(Ex);
                Em *= C1[m+1]   / C1[m];     G = TMP;     Yh.SH(m+2,m+1) += TMP.mxy_matrix(Em);
                Ep *= C2(m+1,m) / C4(m+1,m);              Yh.SH(m+2,m-1) += G.mxy_matrix(Ep);

//              m > 1, 1 < l < l0
                for (size_t l(m+2); l < l0; ++l){
                    MakeGH(l,Yin.SH(l,m));
                    Ex *= A2(l,m) / A1(l-1,m); TMP = H;   Yh.SH(l-1,m  ) += H.mxy_matrix(Ex);
                    Em *= C3[l]   / C1[l-1];   H = TMP;   Yh.SH(l-1,m+1) += TMP.mxy_matrix(Em);
                    Ep *= C4(l,m) / C2(l-1,m);            Yh.SH(l-1,m-1) += H.mxy_matrix(Ep);
                    Ex *= A1(l,m) / A2(l,m);   TMP = G;   Yh.SH(l+1,m  ) += G.mxy_matrix(Ex);
                    Em *= C1[l]   / C3[l];     G = TMP;   Yh.SH(l+1,m+1) += TMP.mxy_matrix(Em);
                    Ep *= C2(l,m) / C4(l,m);              Yh.SH(l+1,m-1) += G.mxy_matrix(Ep); 
                 }
//               m > 1,  l = l0
                 MakeGH(l0,Yin.SH(l0,m));
                 Ex *= A2(l0,m) / A1(l0-1,m); TMP = H;    Yh.SH(l0-1,m  ) += H.mxy_matrix(Ex);
                 Em *= C3[l0]   / C1[l0-1];   H = TMP;    Yh.SH(l0-1,m+1) += TMP.mxy_matrix(Em);
                 Ep *= C4(l0,m) / C2(l0-1,m);             Yh.SH(l0-1,m-1) += H.mxy_matrix(Ep);
             }
         }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         MakeGH(m0,Yin.SH(m0,m0));
         Ep *= C4(m0,m0) / C4(l0,m0-1);              Yh.SH(m0-1,m0-1) += H.mxy_matrix(Ep); 
//       m = m0, l0 > l > m0
         if ( m0 < l0) { 
            Ex *= A1(m0,m0) / A2(l0,m0-1); TMP = G;  Yh.SH(m0+1,m0  ) += TMP.mxy_matrix(Ex); 
            Ep *= C2(m0,m0) / C4(m0,m0);             Yh.SH(m0+1,m0-1) += G.mxy_matrix(Ep);

//          m = m0 , l = m0+1
            MakeGH(m0+1,Yin.SH(m0+1,m0));
            Ex *= A2(m0+1,m0) / A1(m0,m0); TMP = H;        Yh.SH(m0,m0  )  += TMP.mxy_matrix(Ex);
            Ep *= C4(m0+1,m0) / C2(m0,m0);                 Yh.SH(m0,m0-1)  += H.mxy_matrix(Ep);
            if ( m0+1 < l0) { 
                Ex *= A1(m0+1,m0) / A2(m0+1,m0); TMP = G;  Yh.SH(m0+2,m0  )+= TMP.mxy_matrix(Ex);
                Ep *= C2(m0+1,m0) / C4(m0+1,m0);           Yh.SH(m0+2,m0-1)+= G.mxy_matrix(Ep);

//              m = m0, m0+2 < l < l0
                for (size_t l(m0+2); l < l0; ++l){
                    MakeGH(l,Yin.SH(l,m0));
                    Ex *= A2(l,m0) / A1(l-1,m0); TMP = H;  Yh.SH(l-1,m0)   += TMP.mxy_matrix(Ex);
                    Ep *= C4(l,m0) / C2(l-1,m0);           Yh.SH(l-1,m0-1) += H.mxy_matrix(Ep);
                    Ex *= A1(l,m0) / A2(l,  m0); TMP = G;  Yh.SH(l+1,m0  ) += TMP.mxy_matrix(Ex);
                    Ep *= C2(l,m0) / C4(l,  m0);           Yh.SH(l+1,m0-1) += G.mxy_matrix(Ep); 
                 }

//               m > 1,  l = l0
                 MakeGH(l0,Yin.SH(l0,m0));
                 Ex *= A2(l0,m0) / A1(l0-1,m0); TMP = H;    Yh.SH(l0-1,m0)   += TMP.mxy_matrix(Ex);
                 Ep *= C4(l0,m0) / C2(l0-1,m0);             Yh.SH(l0-1,m0-1) += H.mxy_matrix(Ep);
             } 
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void ElectricField::MakeGH(size_t el, SHarmonic& f){
//--------------------------------------------------------------
        GSlice_iter< complex<double> >  it1(f.p0(1));
        complex<double> lc(static_cast< complex<double> >(el)); 
		Axis< complex<double> > invpax(invpr);
		GSlice_iter< complex<double> > iprt(f.p0(el));
		/*if(el == 2){
			G = 0.0;
			H = 0.0;
			return;
		}*/
			
       
		//new version
		/*
        G = f;
		G = G.mpaxis(inv_pls[el-1]);
        G = G.Dp();
		G = G.mpaxis(pls[el-1]);

		H = f;
		H = H.mpaxis(pls[el]);
		H = H.Dp();
		H = H.mpaxis(inv_pls[el]);
		*/

		//Tzou
        invpax *= -2.0*(lc+1.0)*pr.dx();
			
        G = f;                   H = f;
		TMP2 = f;
		H = H.mpaxis(invpax.array());
		invpax *= -lc/(lc+1.0);
		G = G.mpaxis(invpax.array());

        TMP2 = TMP2.Dp();
		H += TMP2;
		G += TMP2;
		

		//mix new G + old H
		/*
        G = f;
		G = G.mpaxis(inv_pls[el-1]);
        G = G.Dp();
		G = G.mpaxis(pls[el-1]);

		H = f;
		TMP2 = f;
		H = H.mpaxis(invpax.array());
        TMP2 = TMP2.Dp();
		H += TMP2;
		*/

		//2nd version
		/*
        G = f;
		G = G.mpaxis(inv_pls[el-1]);
        G = G.Dp();
		G = G.mpaxis(pls[el-1]);

		invpax *= (2.0*lc+1.0)/(lc+1.0);
		H = f;
		H = H.mpaxis(invpax.array());
		H+= G;
*/
		//iprt = G.p0(0);
		//cout<<el<<" "<<*iprt<<endl;
        //GSlice_iter< complex<double> >  itG1(G.p0(1));

        for ( GSlice_iter< complex<double> > itG(G.p0(0)); itG!=itG.end(); ++itG)
            *itG = 0.0;
			//*itG = (*itG1++)*pow(pr(0)/pr(1), el);
		//if(el < 2){
		//iprt = H.p0(4);
		//cout<<el<<" "<<*iprt<<endl;
		size_t ii=0;
			for ( GSlice_iter< complex<double> > itH(H.p0(0)); itH!=itH.end(); ++itH){
				*itH = (*it1++)*Hp0[el];
			}
		//iprt = f.p0(0);
		//cout<<el<<" 0: "<<*iprt<<endl;
		//iprt = f.p0(1);
		//cout<<el<<" 1: "<<*iprt<<endl;
		//iprt = f.p0(2);
		//cout<<el<<" 2: "<<*iprt<<endl;
		//}else{
		//	for ( GSlice_iter< complex<double> > itH(H.p0(0)); itH!=itH.end(); ++itH)
		//		*itH = 0.0;
		//}
		//H=0.0;

		//clean all pi<l / pi<=l
		/*
        for (size_t np=0; np <= el;++np){
		   for ( GSlice_iter< complex<double> > ithl(H.p0(np)); ithl!=ithl.end(); ++ithl){
				   *ithl = 0.0;
		   }
	    }
		*/
		//fits i=l-1, l
		/*
		size_t lp=el+1, lm=el-2, np;
		GSlice_iter< complex<double> > itHH(H.p0(lp));

        for (np=lm; np <lp;++np){
		   //iprt = H.p0(np);
		   itHH = H.p0(lp);
		   if(np>1){
			for ( GSlice_iter< complex<double> > ithl(H.p0(np)); ithl!=ithl.end(); ++ithl){
				    *ithl = (*itHH++)*pow(pr(np)/pr(lp), el-1);
				   // like p0
				   // *ithl = -(*itHH++)*pow(pr(np)/pr(lp), el)/pr(np)*(4.0*lc+2.)*pr.dx();
				//if(el == 4 && abs(*iprt) > 0.0){
				//	cout<<el<<" at "<<ii++<<" t="<<np<<" : "<<(*iprt)<<" from "<<*itHH<<endl;
				//}
				//++iprt;

			}
		   }
	    }
		*/
		//G bound
		/*
		np = el;
		itHH = G.p0(np+1);
		for ( GSlice_iter< complex<double> > ithl(G.p0(np)); ithl!=ithl.end(); ++ithl){
		    // *ithl = (*itHH++)*pow(pr(np)/pr(np+1), el+1);
		    *ithl = 0.0;
		}
*/
		// test
		//H=0.0;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    SHarmonic& ElectricField::G00(SHarmonic& f){ 
//--------------------------------------------------------------
//  Calculation of G00
//--------------------------------------------------------------
        complex<double> p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
                        inv_mp0p1_sq(1.0/(1.0-p0p1_sq)),
                        g_r = -4.0*pr.dx()*pr(0)/(pr(1)*pr(1)),
                        f00,
						g_s = -4.0*pr(0)*pr.dx(),
						inv_p0p1_sq(1.0/(pr(1)*pr(1) - pr(0)*pr(0)));
        GSlice_iter< complex<double> > it0(f.p0(0)), it1(f.p0(1));

        G = f; G = G.Dp();                            

		inv_p0p1_sq *= g_s;
        GSlice_iter< complex<double> >  itG1(G.p0(1));

        for ( GSlice_iter< complex<double> > itG(G.p0(0)); itG!=itG.end(); ++itG){
             f00  =  ( (*it0++) - (*it1)*p0p1_sq)*inv_mp0p1_sq;	
              *itG =  ((*it1++)-f00)*g_r;
			 // simpler
			 //cout<<((*it1) - (*it0))*inv_p0p1_sq/g_s;
			 // *itG = ((*it1++) - (*it0++))*inv_p0p1_sq;
			 //cout<<*itG<<endl;

			 // test gaussian
			 // f00  = (log(*it0)-log(*it1++))/(pr(1)*pr(1) - pr(0)*pr(0));
			 // *itG = (*it0++)*(g_s*f00);

			 // p^l
			 //*itG = (*itG1++) *pr(0)/pr(1);
        }

		//higher order p**2 + p**4
		/*
        GSlice_iter< complex<double> >  it2(f.p0(2));
		complex<double> b,c, gs = -4.0*pr.dx();
		complex<double> p02=pow(pr(0),2),p04=pow(pr(0),4);
		complex<double> p12=pow(pr(1),2),p14=pow(pr(1),4);
		complex<double> p22=pow(pr(2),2),p24=pow(pr(2),4);
		complex<double> br= (p12-p02)*(p24-p04) - (p22-p02)*(p14-p04);
		complex<double> cr= (p22-p02)*(p14-p04) - (p12-p02)*(p24-p04);

        for ( GSlice_iter< complex<double> > itG(G.p0(0)); itG!=itG.end(); ++itG){
			b = ((*it1 - *it0)*(p24-p04) - (*it2 - *it0)*(p14-p04))/br;
			c = ((*it1 - *it0)*(p22-p02) - (*it2 - *it0)*(p12-p02))/cr;
			*itG = (b*pr(0) + 2.0*c*pow(pr(0),3))*gs;
			//cout<<b<<" and  "<<c<<" gives "<<*itG<<endl;
			++it0;
			++it1;
			++it2;
		}
		*/
        return G;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Magnetic Field
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    MagneticField::MagneticField(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(m0+1),      B1(l0+1),
         A2(l0+1,m0+1),
         A3(0.5),
         FLM(Yslope.SH(0,0)) 
         {
         complex<double> lc, mc;
         complex<double> c01(0.0,1.0);

//       Calculate the "A1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t m=0; m<m0+1; ++m){
             mc = static_cast< complex<double> >(m);
             A1[m] = (-1.0)*c01*mc;
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l)
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 A2(l,m) = (-0.5)*(lc+1.0-mc)*(lc+mc);
             }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         A3 = 0.5;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
            lc = static_cast< complex<double> >(l);
            B1[l] = (-1.0)*lc*(lc+1.0); 
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& MagneticField::Bxyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        //Matrix2D< complex<double> > Bx(Yin.EMF().Bx().matrix());
        //Matrix2D< complex<double> > Bm(Yin.EMF().By().matrix()); Bm *= (-1.0)*ii; Bm += Yin.EMF().Bz().matrix();
        //Matrix2D< complex<double> > Bp(Yin.EMF().By().matrix()); Bp *=  ii; Bp += Yin.EMF().Bz().matrix();
		// --- e-i
        Matrix2D< complex<double> > Bx(Yin.EMF().Bx().matrix()); Bx *= (-1.0);
        Matrix2D< complex<double> > Bm(Yin.EMF().By().matrix()); Bm *= ii;        Bm -= Yin.EMF().Bz().matrix();
        Matrix2D< complex<double> > Bp(Yin.EMF().By().matrix()); Bp *= (-1.0)*ii; Bp -= Yin.EMF().Bz().matrix();

//      m = 0, 1 <= l < l0+1                                    //WHY? m>0 in P15
        Bp *= A3;
        for (size_t l(1); l < l0+1; ++l){
            FLM = Yin.SH(l,0);      Yh.SH(l,1) += FLM.mxy_matrix(Bp);
        }

//      m = 1, l = 1
        FLM = Yin.SH(1,1); Bx *= A1[1];                           Yh.SH(1,1) += FLM.mxy_matrix(Bx); 
        FLM = Yin.SH(1,1); Bm *= B1[1]; FLM = FLM.mxy_matrix(Bm); Yh.SH(1,0) += FLM.Re(); 
//      m = 1, l > 1
        for (size_t l(2); l < l0+1; ++l){
            FLM = Yin.SH(l,1);                                                Yh.SH(l,2) += FLM.mxy_matrix(Bp);     //WHY- Bp,Bx,Bm already *A1, A3,B1,etc above 
            FLM = Yin.SH(l,1);                                                Yh.SH(l,1) += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(l,1); Bm *= B1[l]/B1[l-1]; FLM = FLM.mxy_matrix(Bm); Yh.SH(l,0) += FLM.Re(); 
        }        
        Bm *= 1.0/B1[l0];        

//      m > 1, l = m
        for (size_t m(2); m < m0; ++m){
            FLM = Yin.SH(m,m); Bx *= A1[m]/A1[m-1];                   Yh.SH(m,m  ) += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(m,m); Bm *= A2(m,m);                         Yh.SH(m,m-1) += FLM.mxy_matrix(Bm); 
            for (size_t l(m+1); l < l0+1; ++l){
                FLM = Yin.SH(l,m);                                    Yh.SH(l,m+1) += FLM.mxy_matrix(Bp); 
                FLM = Yin.SH(l,m);                                    Yh.SH(l,m  ) += FLM.mxy_matrix(Bx); 
                FLM = Yin.SH(l,m); Bm *= A2(l,m)/A2(l-1,m);           Yh.SH(l,m-1) += FLM.mxy_matrix(Bm); 
            }
            Bm *= 1.0/A2(l0,m);
        }

//      m = m0, l >= m0
        FLM = Yin.SH(m0,m0); Bx *= A1[m0]/A1[m0-1];                   Yh.SH(m0,m0)   += FLM.mxy_matrix(Bx); 
        FLM = Yin.SH(m0,m0); Bm *= A2(m0,m0)/*/A2(l0,m0-1)*/;             Yh.SH(m0,m0-1) += FLM.mxy_matrix(Bm);
        for (size_t l(m0+1); l < l0+1; ++l){
            FLM = Yin.SH(l,m0);                                     Yh.SH(l,m0  )  += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(l,m0); Bm *= A2(l,m0)/A2(l-1,m0);          Yh.SH(l,m0-1)  += FLM.mxy_matrix(Bm); 
        }

        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Current in the x direction
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Current::Current(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         jayx(Yslope.FLD(0)), 
         jayy(Yslope.FLD(0)), 
         jayz(Yslope.FLD(0)), 
         p3og(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         Delta_p(static_cast< complex<double> >(Inputdata::IN().inp().pr.dx())),
         small(static_cast< complex<double> >(Inputdata::IN().inp().pr(0)))
         {
             for (size_t i(0); i<p3og.dim(); ++i) // p^3/gamma 
                p3og(i) = (p3og(i)*p3og(i))*(p3og(i)/sqrt(1.0+p3og(i)*p3og(i))); 

             small *= small; small *= small; small *= Inputdata::IN().inp().pr(0);      //WHY- P12 for [0,p0] 
             small *= 0.2; small *= 1.0/(Inputdata::IN().inp().pr(1)); 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Current::Jxyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> c01(0.0,1.0);
        SHarmonic f10(Yin.SH(1,0)); 
        SHarmonic f11(Yin.SH(1,1)); 
        f10  = f10.mpaxis(p3og.array());
        f11  = f11.mpaxis(p3og.array());

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      Basic integration
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy)  = f10(0,ix,iy);
                jayx(ix,iy) += f10(f10.nump()-1,ix,iy);
                jayx(ix,iy) *= 0.5;                                     //WHY? just *0.5 for p1 & pmax?
                for (size_t ip (1); ip < f10.nump()-1; ++ip) jayx(ix,iy)+= f10(ip,ix,iy);
            }
        }
        jayx *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy) += small*f10(1,ix,iy);
             }
        }

//      Calculate the final minus current
        jayx *= 4.0 * M_PI / 3.0; 

        Yh.EMF().Ex() += jayx;// (minus minus ---> +current)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//      Basic integration
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy)  = f11(0,ix,iy);
                jayz(ix,iy) += f11(f11.nump()-1,ix,iy);
                jayz(ix,iy) *= 0.5;
                for (size_t ip (1); ip < f11.nump()-1; ++ip) jayz(ix,iy)+= f11(ip,ix,iy);
            }
        }
        jayz *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy) += small*f11(1,ix,iy);
             }
        }

//      Calculate the final minus complex current 
        jayz *= 8.0 * M_PI / 3.0;

        jayy = jayz;      jayy  = jayy.Re();
        jayz -= jayy;     jayz *= c01;
        
//      minus minus ---> +current        
        Yh.EMF().Ey() += jayy;
        Yh.EMF().Ez() += jayz;

        return Yh;
        
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Spatial Advection 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Spatial_Advection::Spatial_Advection(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         fd1(Yslope.SH(0,0)), fd2(Yslope.SH(0,0)),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(l0+1,m0+1), A2(l0+1,m0+1),
         B1(l0+1),      B2(l0+1),
         C1(l0+1),      C3(l0+1),
         C2(l0+1,m0+1), C4(l0+1,m0+1),
         vr(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx()))
         {

         for (size_t i(0); i<vr.dim(); ++i) vr(i) = vr(i)/(sqrt(1.0+vr(i)*vr(i)));  //WHY- p=gamma*m*v
         idx = -1.0/(2.0*idx);                                  //WHY- using central diff 
         idy = -1.0/(2.0*idy); 

         complex<double> lc, mc;

//       Calculate the "A1, A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 A1(l,m) = idx *(-1.0) * (lc-mc+1.0) / (2.0*lc+1.0);
                 A2(l,m) = idx *(-1.0) * (lc+mc)     / (2.0*lc+1.0);    
             }
         }

//         A1(0,0) = -1.0;
//         A1 *= 0.5 / pr.dx(); A2 *= 0.5 / pr.dx();
//         A2(0,0) = 1.0;
//         for (size_t m(1); m < m0+1; ++m) A2(m,m) = A2(l0,m-1);
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1, B2" parameters                      //WHY? not just real
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             B1[l] = idy * (lc + 1.0) * lc / (2.0*lc + 1.0);
             B2[l] = (-1.0)*B1[l];
         }

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C1, C3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             C1[l] = (-0.5) * idy / (2.0*lc + 1.0);
             C3[l] = (-1.0) * C1[l];
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C2, C4" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 C2(l,m) = idy * 0.5 * (lc + 2.0 - mc)*(lc - mc + 1.0) / (2.0*lc + 1.0);
                 C4(l,m) = idy * (-0.5) * (lc + mc - 1.0)*(lc + mc) / (2.0*lc + 1.0);
             }
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Spatial_Advection::Axy(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the spatial advection
//--------------------------------------------------------------
        Axis< complex<double> > vt(vr);        

//      Advection in x
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t m(0); m < ((m0<l0)?(m0+1):m0); ++m){
            fd1 = Yin.SH(m,m);     fd1 = fd1.Dx();
            vt *= A1(m,m);                          Yh.SH(m+1,m) += fd1.mpaxis(vt.array());
            for (size_t l(m+1); l < l0; ++l) {
               fd1 = Yin.SH(l,m);          fd1 = fd1.Dx();  
               vt *= A2(l,m)/A1(l-1,m);    fd2 = fd1;  Yh.SH(l-1,m) += fd1.mpaxis(vt.array());  
               vt *= A1(l,m)/A2(l  ,m);                Yh.SH(l+1,m) += fd2.mpaxis(vt.array());  
            }
            fd1 = Yin.SH(l0,m);            fd1 = fd1.Dx();
            vt *= A2(l0,m)/A1(l0-1,m);                 Yh.SH(l0-1,m) += fd1.mpaxis(vt.array());  
            vt *= 1.0/A2(l0,m); 
        }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
            
//       m = 0, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(0,0);     fd1 = fd1.Dy();                                                 //WHY? when it has df00/dy?
         vt *= C1[0];                            Yh.SH(1,1) += fd1.mpaxis(vt.array());
         fd1 = Yin.SH(1,0);     fd1 = fd1.Dy();
         vt *= C1[1]/C1[0];                      Yh.SH(2,1) += fd1.mpaxis(vt.array());
         for (size_t l(2); l < l0; ++l) {
             fd1 = Yin.SH(l,0);     fd1 = fd1.Dy();
             vt *= C3[l]/C1[l-1];   fd2 = fd1;   Yh.SH(l-1,1) += fd1.mpaxis(vt.array());
             vt *= C1[l]/C3[l];                  Yh.SH(l+1,1) += fd2.mpaxis(vt.array());        //WHY? just *-1
         }
         fd1 = Yin.SH(l0,0);    fd1 = fd1.Dy();
         vt *= C3[l0]/C1[l0-1];              Yh.SH(l0-1,1) += fd1.mpaxis(vt.array());
         vt *= 1.0   /C3[l0];
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

         
//       m = 1, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(1,1);     fd1 = fd1.Dy();
         vt *= C1[1];           fd2 = fd1;                                 Yh.SH(2,2) += fd1.mpaxis(vt.array());
         vt *= B2[1]/C1[1];     fd1 = fd2;  fd2 = fd2.mpaxis(vt.array());  Yh.SH(0,0) += fd2.Re();
         vt *= B1[1]/B2[1];                 fd1 = fd1.mpaxis(vt.array());  Yh.SH(2,0) += fd1.Re();

         fd1 = Yin.SH(2,1);     fd1 = fd1.Dy();
         vt *= C1[2]/B1[1];     fd2 = fd1;                                 Yh.SH(3,2) += fd1.mpaxis(vt.array());
         vt *= B2[2]/C1[2];     fd1 = fd2;  fd2 = fd2.mpaxis(vt.array());  Yh.SH(1,0) += fd2.Re();
         vt *= B1[2]/B2[2];                 fd1 = fd1.mpaxis(vt.array());  Yh.SH(3,0) += fd1.Re();
       
         for (size_t l(3); l < l0; ++l){                                                                            //WHY- coz C3:l-1,m+1
             fd1 = Yin.SH(l,1);     fd1 = fd1.Dy();                           
             vt *= C3[l]/B1[l-1];   fd2 = fd1;                                Yh.SH(l-1,2) += fd1.mpaxis(vt.array());
             vt *= C1[l]/C3[l];     fd1 = fd2;                                Yh.SH(l+1,2) += fd2.mpaxis(vt.array());
             vt *= B2[l]/C1[l];     fd2 = fd1;  fd1 = fd1.mpaxis(vt.array()); Yh.SH(l-1,0) += fd1.Re();
             vt *= B1[l]/B2[l];                 fd2 = fd2.mpaxis(vt.array()); Yh.SH(l+1,0) += fd2.Re();
         }
         fd1 = Yin.SH(l0,1);     fd1 = fd1.Dy();                           
         vt *= C3[l0]/B1[l0-1];  fd2 = fd1;                                Yh.SH(l0-1,2) += fd1.mpaxis(vt.array());
         vt *= B2[l0]/C3[l0];                fd2 = fd2.mpaxis(vt.array()); Yh.SH(l0-1,0) += fd2.Re();
         vt *= 1.0   /B2[l0];
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
//       m > 1, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t m(2); m < m0; ++m){
//           m > 1, l = m
             fd1 = Yin.SH(m,m);       fd1 = fd1.Dy();
             vt *= C4(m,m);           fd2 = fd1;     Yh.SH(m-1,m-1) += fd1.mpaxis(vt.array());
             vt *= C2(m,m)/C4(m,m);   fd1 = fd2;     Yh.SH(m+1,m-1) += fd2.mpaxis(vt.array());
             vt *= C1[m]  /C2(m,m);                  Yh.SH(m+1,m+1) += fd1.mpaxis(vt.array());

//           m > 1, l = m+1
             fd1 = Yin.SH(m+1,m);       fd1 = fd1.Dy();
             vt *= C4(m+1,m)/C1[m];     fd2 = fd1;     Yh.SH(m  ,m-1) += fd1.mpaxis(vt.array());
             if (m+1 < l0) {
                 vt *= C2(m+1,m)/C4(m+1,m);   fd1 = fd2;   Yh.SH(m+2,m-1) += fd2.mpaxis(vt.array());
                 vt *= C1[m+1]  /C2(m+1,m);                Yh.SH(m+2,m+1) += fd1.mpaxis(vt.array());

//               m > 1, 3 < l < l0
                 for (size_t l(m+2); l < l0; ++l){
                     fd1 = Yin.SH(l,m);       fd1 = fd1.Dy();
                     vt *= C4(l,m)/C1[m+1];   fd2 = fd1;     Yh.SH(l-1,m-1) += fd1.mpaxis(vt.array());
                     vt *= C2(l,m)/C4(l,m);   fd1 = fd2;     Yh.SH(l+1,m-1) += fd2.mpaxis(vt.array());
                     vt *= C3[l]  /C2(l,m);   fd2 = fd1;     Yh.SH(l-1,m+1) += fd1.mpaxis(vt.array());
                     vt *= C1[l]  /C3[l];                    Yh.SH(l+1,m+1) += fd2.mpaxis(vt.array());
                 }
                 
                 fd1 = Yin.SH(l0,m);       fd1 = fd1.Dy();
                 vt *= C4(l0,m)/C1[l0-1];  fd2 = fd1;       Yh.SH(l0-1,m-1) += fd1.mpaxis(vt.array());
                 vt *= C3[l0]/C4(l0,m);                     Yh.SH(l0-1,m+1) += fd2.mpaxis(vt.array());
                 vt *= 1.0/C3[l0];
             }
             else {
                 vt *= 1.0/C4(m+1,m);
             }
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(m0,m0);       fd1 = fd1.Dy();
         vt *= C4(m0,m0);           fd2 = fd1;     Yh.SH(m0-1,m0-1) += fd1.mpaxis(vt.array());
         if (m0 < l0) {
             vt *= C2(m0,m0)/C4(m0,m0);            Yh.SH(m0+1,m0-1) += fd2.mpaxis(vt.array());
             for (size_t l(m0+1); l < l0; ++l){
                 fd1 = Yin.SH(l,m0);          fd1 = fd1.Dy();
                 vt *= C4(l,m0)/C2(l-1,m0);   fd2 = fd1;     Yh.SH(l-1,m0-1) += fd1.mpaxis(vt.array());
                 vt *= C2(l,m0)/C4(l,  m0);                  Yh.SH(l+1,m0-1) += fd2.mpaxis(vt.array());
             }
             fd1 = Yin.SH(l0,m0);           fd1 = fd1.Dy();
             vt *= C4(l0,m0)/C2(l0-1,m0);   Yh.SH(l0-1,m0-1) += fd1.mpaxis(vt.array());
         }

        return Yh;    

    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for Maxwell's Equations
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    MaxwellEq::MaxwellEq(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         tmpEB(Yslope.FLD(0)), 
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx()))
         {
         idx = -1.0/(2.0*idx); 
         idy = -1.0/(2.0*idy); 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& MaxwellEq::Faraday(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for Faraday's Law 
//--------------------------------------------------------------

//      dBx/dt += - dEz/dy  
        tmpEB          = Yin.EMF().Ez(); 
        tmpEB         *= (-1.0) * idy;
        Yh.EMF().Bx() += tmpEB.Dy();

//      dBy/dt +=   dEz/dx       
        tmpEB          = Yin.EMF().Ez(); 
        tmpEB         *= idx;
        Yh.EMF().By() += tmpEB.Dx();        

//      dBz/dt +=   dEx/dy       
        tmpEB          = Yin.EMF().Ex(); 
        tmpEB         *= idy;
        Yh.EMF().Bz() += tmpEB.Dy();    

//      dBz/dt += - dEy/dx       
        tmpEB          = Yin.EMF().Ey(); 
        tmpEB         *= (-1.0) * idx;
        Yh.EMF().Bz() += tmpEB.Dx();   

        return Yh;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Stat& MaxwellEq::Ampere(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for Ampere's Law 
//--------------------------------------------------------------

//      dEx/dt +=   dBz/dy       
        tmpEB          = Yin.EMF().Bz(); 
        tmpEB         *= idy;
        Yh.EMF().Ex() += tmpEB.Dy();

//      dEy/dt +=  - dBz/dx       
        tmpEB          = Yin.EMF().Bz(); 
        tmpEB         *= (-1.0) * idx;
        Yh.EMF().Ey() += tmpEB.Dx();        

//      dEz/dt +=  - dBx/dy       
        tmpEB          = Yin.EMF().Bx(); 
        tmpEB         *= (-1.0) * idy;
        Yh.EMF().Ez() += tmpEB.Dy();    

//      dEz/dt += dBy/dx       
        tmpEB          = Yin.EMF().By(); 
        tmpEB         *=  idx;
        Yh.EMF().Ez() += tmpEB.Dx();   

        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for Ohm
//**************************************************************
//**************************************************************
//
//**************************************************************
//--------------------------------------------------------------
	OHMic::OHMic(Stat& Yin)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yin),
         tmpEB(Yin.FLD(0)), 
		 // --- e-i
         Nbc(Inputdata::IN().list().RKLevel),
		 pr(Inputdata::IN().inp().pr),
		 gamma(pr),
		 invgamma(pr),
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx())),
		 // for AD use, integrals
         vr(0.0, Inputdata::IN().list().nump), 
         U4(0.0, Inputdata::IN().list().nump), 
         U4m1(0.0, Inputdata::IN().list().nump), 
         U2(0.0, Inputdata::IN().list().nump), 
         U2m1(0.0, Inputdata::IN().list().nump)
         {
         idx = -1.0/(2.0*idx); 
         idy = -1.0/(2.0*idy); 
		 //--- e-i version
		 p0p1_sq      = pr(0)*pr(0)/(pr(1)*pr(1));
		 inv_mp0p1_sq = 1.0/(1.0-p0p1_sq);
         for (size_t ip(0); ip < pr.dim(); ++ip) gamma(ip)    = sqrt(1.0+pr(ip)*pr(ip));
         for (size_t ip(0); ip < pr.dim(); ++ip) {
			 invgamma(ip) = 1.0 / gamma(ip);
             //vr[ip] = pr(ip) * invgamma(ip);
		 }

		 // for central difference, needs +2
         szx = Inputdata::IN().inp().x.dim()-2*Nbc+2;
         szy = Inputdata::IN().inp().y.dim()-2*Nbc+2;
		 //--- e-i
		 /*double tmp = 2.0*M_PI;//4pi*0.5
         for (size_t i(1); i < U4.size(); ++i) {
             U4[i]   = tmp * pow(vr[i],4)     * (vr[i]-vr[i-1]);         
             U4m1[i] = tmp * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);         
             U2[i]   = tmp * vr[i]*vr[i]      * (vr[i]-vr[i-1]);         
             U2m1[i] = tmp * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);         
         }
		 */
    }
//--------------------------------------------------------------
//	Ampere for J
//--------------------------------------------------------------
    void OHMic::J_Ampere(Stat& Yi, Stat& Ye){ 
//      Jx =   dBz/dy       
        tmpEB          = Yi.EMF().Bz(); 
        tmpEB         *= idy;
        Ye.EMF().Ex()  = tmpEB.Dy();

//      Jy =  - dBz/dx       
        tmpEB          = Yi.EMF().Bz(); 
        tmpEB         *= (-1.0) * idx;
        Ye.EMF().Ey()  = tmpEB.Dx();

//      Jz =  - dBx/dy       
        tmpEB          = Yi.EMF().Bx(); 
        tmpEB         *= (-1.0) * idy;
        Ye.EMF().Ez()  = tmpEB.Dy();    

//      Jz += dBy/dx       
        tmpEB          = Yi.EMF().By(); 
        tmpEB         *=  idx;
        Ye.EMF().Ez() += tmpEB.Dx();   

        return;
    }
//--------------------------------------------------------------
//	for density
//--------------------------------------------------------------
    void OHMic::J_density(Stat& Yi, Stat& Ye, const int st){ 
		//st = 0: initial, st=1: in RK steps

		SHarmonic& F0(Yi.SH(0,0));
		EMF2D& VD(Ye.EMF());
        Matrix2D<double> density(szx,szy);  
		complex<double> f00;
		int yid, xid;

		//N_i
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;

                density(ix,iy)  = F0(0,xid,yid).real()*pr(0)*pr(0);
                density(ix,iy) += F0(F0.nump()-1,xid,yid).real()*pr(pr.dim()-1)*pr(pr.dim()-1);
                density(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F0.nump()-1; ++ip){ 
                    density(ix,iy) += pr(ip)*pr(ip)*F0(ip,xid,yid).real();
                }

                f00  =  F0(1,xid,yid)/5.0*p0p1_sq;
                f00 += (F0(0,xid,yid) - F0(1,xid,yid)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/3.0-1.0/5.0*p0p1_sq); 
                //f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 
                f00 *= pow(pr(0),3)/pr.dx();

                density(ix,iy) += f00.real();
				density(ix,iy) *= 4.0 * M_PI * pr.dx();

                if (density(ix,iy) < (4.0*FLT_EPSILON))
					VD.By()(xid,yid).real() = 1.0;
				else
					VD.By()(xid,yid).real() = density(ix,iy);

            }
        }
		//N_e
        //density *= pr.dx()*Inputdata::IN().list().Zeta;

        return;
    }
//--------------------------------------------------------------
//	Ve = (Vi-J/ZN_i), J/n is the Hall
//	E  = -VxB
//--------------------------------------------------------------
    void OHMic::J_E(Stat& Yi, Stat& Ye){ 
        // ---------  Vx + Vyz  ------------------------------------ 
		SHarmonic& F10(Yi.SH(1,0)),F11(Yi.SH(1,1));
		EMF2D& EF(Yi.EMF()), VD(Ye.EMF());
        Matrix2D<double> Vx(szx,szy);  
        Matrix2D<complex<double> > Vyz(szx,szy); 
		complex<double> f00, f11;
		double vx, vy, vz, tmp;;
		int yid, xid;

		// p
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
				//px
                Vx(ix,iy)  = F10(0,xid,yid).real() * pow(pr(0),3)*invgamma(0);
                Vx(ix,iy) += F10(F10.nump()-1,xid,yid).real() *
                                  pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                Vx(ix,iy) *= 0.5;
				//pyz
                Vyz(ix,iy)  = F11(0,xid,yid) * pow(pr(0),3)*invgamma(0);
                Vyz(ix,iy) += F11(F11.nump()-1,xid,yid) * 
                                 pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                Vyz(ix,iy) *= 0.5;
				//ingetral
                for (size_t ip(1); ip < F10.nump()-1; ++ip) { 
                   Vx(ix,iy) += pow(pr(ip),3)*invgamma(ip)*F10(ip,xid,yid).real();

                   Vyz(ix,iy)+= pow(pr(ip),3)*invgamma(ip)*F11(ip,xid,yid);
                }
				//0->p0
                f00  = F10(1,xid,yid);
                f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                Vx(ix,iy) += f00.real();

                f11  =  F11(1,xid,yid);
                f11 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                Vyz(ix,iy) += f11;

            }
        }
        Vx  *= (4.0 * M_PI / 3.0) * pr.dx();
        Vyz *= (8.0 * M_PI / 3.0) * pr.dx();

		//pe
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
				// Divide by the density to get ve = vi - J/n = Vi/Ni - J/Ni
				tmp= VD.By()(xid,yid).real();
				// MHD .or.
                vx = Vx(ix,iy)        /tmp;
                vy = Vyz(ix,iy).real()/tmp;
                vz =-Vyz(ix,iy).imag()/tmp;
				// .or. Hall
                //vx = Vx(ix,iy)        /tmp - VD.Ex()(xid,yid).real();
                //vy = Vyz(ix,iy).real()/tmp - VD.Ey()(xid,yid).real();
                //vz =-Vyz(ix,iy).imag()/tmp - VD.Ez()(xid,yid).real();

				// J/ne
                VD.Ex()(xid,yid) /= (Inputdata::IN().list().Zeta*tmp);
                VD.Ey()(xid,yid) /= (Inputdata::IN().list().Zeta*tmp);
                VD.Ez()(xid,yid) /= (Inputdata::IN().list().Zeta*tmp);

				// -Vi x B
                EF.Ex()(xid,yid).real()  = (vz*EF.By()(xid,yid).real()-vy*EF.Bz()(xid,yid).real());
                EF.Ey()(xid,yid).real()  = (vx*EF.Bz()(xid,yid).real()-vz*EF.Bx()(xid,yid).real());
                EF.Ez()(xid,yid).real()  = (vy*EF.Bx()(xid,yid).real()-vx*EF.By()(xid,yid).real());

				// -Vi x B - grad.Pe
                //EF.Ex()(xid,yid).real()  = (vz*EF.By()(xid,yid).real()-vy*EF.Bz()(xid,yid).real()) + VD.Ex()(xid,yid).imag();
                //EF.Ey()(xid,yid).real()  = (vx*EF.Bz()(xid,yid).real()-vz*EF.Bx()(xid,yid).real()) + VD.Ey()(xid,yid).imag();
                //EF.Ez()(xid,yid).real()  = (vy*EF.Bx()(xid,yid).real()-vx*EF.By()(xid,yid).real());

            }
        }

        return;
    }
//--------------------------------------------------------------
// J_Ampere + J_density, <|v|^2>, not <|v-V|^2>
//--------------------------------------------------------------
    void OHMic::J_AD(Stat& Yi, Stat& Ye){ 
		SHarmonic& F0(Yi.SH(0,0));
		EMF2D& VD(Ye.EMF());
        Matrix2D<double> density(szx,szy);
		Field2D Pp(Inputdata::IN().inp().x.dim(),Inputdata::IN().inp().y.dim()); 
		complex<double> f00;
		int yid, xid;
		double tmp_den, tmp_pe, fn, fn1;

		// Ne, Pe
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;

				tmp_den = 0.0;
				tmp_pe  = 0.0;
                for (size_t ip(1); ip < F0.nump(); ++ip) {
					fn  = F0(ip,xid,yid).real();
					fn1 = F0(ip-1,xid,yid).real();
					tmp_den += U2[ip]*fn+U2m1[ip]*fn1;
					tmp_pe  += U4[ip]*fn+U4m1[ip]*fn1;
                }
                density(ix,iy) = tmp_den;	   //Ni
				// vt**2/3?
				Pp(xid,yid)      = tmp_pe; //Ni*VTi**2
				//cout<<xid<<" "<<yid<<" Pp:"<<tmp_pe<<endl;
            }
        }

		// -grad.Pp
		Pp      *= (-1.0)*Inputdata::IN().list().Zeta;
		tmpEB    = Pp;
        tmpEB   *= idx;
        VD.Bx()  = tmpEB.Dx();
		tmpEB    = Pp;
        tmpEB   *= idy;
        VD.Bz()  = tmpEB.Dy();

		// J_Ampere
		// Jx =   dBz/dy       
        tmpEB    = Yi.EMF().Bz(); 
        tmpEB   *= idy;
        VD.Ex()  = tmpEB.Dy();

		// Jy =  - dBz/dx       
        tmpEB    = Yi.EMF().Bz(); 
        tmpEB   *= (-1.0) * idx;
        VD.Ey()  = tmpEB.Dx();

		// Jz =  - dBx/dy       
        tmpEB    = Yi.EMF().Bx(); 
        tmpEB   *= (-1.0) * idy;
        VD.Ez()  = tmpEB.Dy();    

		// Jz +=   dBy/dx       
        tmpEB    = Yi.EMF().By(); 
        tmpEB   *= idx;
        VD.Ez() += tmpEB.Dx();   


        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
				// Ni
				tmp_den = density(ix,iy);
				VD.By()(xid,yid).real() = tmp_den;
				// grad.Pe / (Mi.Ni)
				//VD.Ex()(xid,yid).imag() = VD.Bx()(xid,yid).real()/density(ix,iy);
				//VD.Ey()(xid,yid).imag() = VD.Bz()(xid,yid).real()/density(ix,iy);
				VD.Ex()(xid,yid).imag() = VD.Bx()(xid,yid).real();
				VD.Ey()(xid,yid).imag() = VD.Bz()(xid,yid).real();
				// J/Ni
                //VD.Ex()(xid,yid).real() /= tmp;
                //VD.Ey()(xid,yid).real() /= tmp;
                //VD.Ez()(xid,yid).real() /= tmp;
				
                VD.Ex()(xid,yid)        /= tmp_den;
                VD.Ey()(xid,yid)        /= tmp_den;
                VD.Ez()(xid,yid).real() /= tmp_den;
            }
        }

        return;
    }
//--------------------------------------------------------------
// Tzoufras' pressure term, Pe = Pi
//--------------------------------------------------------------
    void OHMic::J_Pre(Stat& Yi, Stat& Ye){ 

		SHarmonic& F0(Yi.SH(0,0));
		EMF2D& VD(Ye.EMF());
        Matrix2D<double> density(szx,szy);
		complex<double> f00;
		int yid, xid;

		// Ni, J_density ---------------------------------------------------
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;

                density(ix,iy)  = F0(0,xid,yid).real()*pr(0)*pr(0);
                density(ix,iy) += F0(F0.nump()-1,xid,yid).real()*pr(pr.dim()-1)*pr(pr.dim()-1);
                density(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F0.nump()-1; ++ip){ 
                    density(ix,iy) += pr(ip)*pr(ip)*F0(ip,xid,yid).real();
                }

                f00  =  F0(1,xid,yid)/5.0*p0p1_sq;
                f00 += (F0(0,xid,yid) - F0(1,xid,yid)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/3.0-1.0/5.0*p0p1_sq); 
                //f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 
                f00 *= pow(pr(0),3)/pr.dx();

                density(ix,iy) += f00.real();
				density(ix,iy) *= 4.0 * M_PI * pr.dx();

                if (density(ix,iy) < (4.0*FLT_EPSILON))
					VD.By()(xid,yid).real() = 1.0;
				else
					VD.By()(xid,yid).real() = density(ix,iy);

            }
        }

        // --------- J_Ampere -----------------------------------------
//      Jx =   dBz/dy       
        tmpEB          = Yi.EMF().Bz(); 
        tmpEB         *= idy;
        Ye.EMF().Ex()  = tmpEB.Dy();

//      Jy =  - dBz/dx       
        tmpEB          = Yi.EMF().Bz(); 
        tmpEB         *= (-1.0) * idx;
        Ye.EMF().Ey()  = tmpEB.Dx();

//      Jz =  - dBx/dy       
        tmpEB          = Yi.EMF().Bx(); 
        tmpEB         *= (-1.0) * idy;
        Ye.EMF().Ez()  = tmpEB.Dy();    

//      Jz += dBy/dx       
        tmpEB          = Yi.EMF().By(); 
        tmpEB         *=  idx;
        Ye.EMF().Ez() += tmpEB.Dx();   

		// Vi, J_E ---------------------------------------------------------
		SHarmonic& F10(Yi.SH(1,0)),F11(Yi.SH(1,1));
		EMF2D& EF(Yi.EMF());
        Matrix2D<double> Vx(szx,szy);  
        Matrix2D<complex<double> > Vyz(szx,szy); 
		complex<double> f11;

        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
				//px
                Vx(ix,iy)  = F10(0,xid,yid).real() * pow(pr(0),3)*invgamma(0);
                Vx(ix,iy) += F10(F10.nump()-1,xid,yid).real() *
                                  pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                Vx(ix,iy) *= 0.5;
				//pyz
                Vyz(ix,iy)  = F11(0,xid,yid) * pow(pr(0),3)*invgamma(0);
                Vyz(ix,iy) += F11(F11.nump()-1,xid,yid) * 
                                 pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                Vyz(ix,iy) *= 0.5;
				//integrate
                for (size_t ip(1); ip < F10.nump()-1; ++ip) { 
                   Vx(ix,iy) += pow(pr(ip),3)*invgamma(ip)*F10(ip,xid,yid).real();

                   Vyz(ix,iy)+= pow(pr(ip),3)*invgamma(ip)*F11(ip,xid,yid);
                }
				//0->p0
                f00  = F10(1,xid,yid);
                f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                Vx(ix,iy) += f00.real();

                f11  =  F11(1,xid,yid);
                f11 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                Vyz(ix,iy) += f11;

            }
        }
        Vx  *= (4.0 * M_PI / 3.0) * pr.dx();
        Vyz *= (8.0 * M_PI / 3.0) * pr.dx();

        // --------- VxVx+VyVy+VzVz -----------------------------------	  //WHY- Note(137)
		// WHY- Note(137)+(140)+(143) = 4pi*f00
		// now density = moment_1

        // integral for "f_0^0"
        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;

                density(ix,iy)  = F0(0,xid,yid).real() * pow(pr(0),4) *
                                      pow(invgamma(0),2);
                density(ix,iy) += F0(F0.nump()-1,xid,yid).real() * pow(pr(pr.dim()-1),4)* 
                                      pow(invgamma(pr.dim()-1),2);
                density(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F0.nump()-1; ++ip){ 
					density(ix,iy) += pow(pr(ip),4)*pow(invgamma(ip),2)*
                                      F0(ip,xid,yid).real();
                }

                f00  =  F0(1,xid,yid)/7.0*p0p1_sq;
                f00 += (F0(0,xid,yid) - F0(1,xid,yid)*p0p1_sq)*inv_mp0p1_sq
                      *(1.0/5.0-1.0/7.0*p0p1_sq); 
                f00 *= pow(pr(0),4)*(pr(0)/pr.dx()); 
                density(ix,iy) += f00.real();

            }
        }
        density *= 4.0 * M_PI * pr.dx();


        // ---------  v^2 -----------------------------------------
		// save data to Ye
		Field2D Pp(Inputdata::IN().inp().x.dim(),Inputdata::IN().inp().y.dim()); Pp=0.0; 
        Matrix2D<double> Vsq(szx,szy);

        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
	
				// VxVx + VyVy + VzVz - V^2
                Pp(xid,yid).real() = density(ix,iy)
					               -(Vx(ix,iy)*Vx(ix,iy)
									+Vyz(ix,iy).real()*Vyz(ix,iy).real()
									+Vyz(ix,iy).imag()*Vyz(ix,iy).imag() );
				Vsq(ix,iy)         = Pp(xid,yid).real()/VD.By()(xid,yid).real();
				//--- for adiab
                Pp(xid,yid).real() = pow(VD.By()(xid,yid).real(), 5.0/3.0)*Ye.SH(0,0)(0,xid,yid).real();
                //Pp(xid,yid).real() = Ye.SH(0,0)(0,xid,yid).real();

            }
        }

        // ---------  -grad.Pre -----------------------------------------
		//Pp      *= (-1.0/3.0);

		//--- for adiab
		Pp      *= (-1.0);
		
		tmpEB    = Pp;
        tmpEB   *= idx;
        VD.Bx()  = tmpEB.Dx();
		tmpEB    = Pp;
        tmpEB   *= idy;
        VD.Bz()  = tmpEB.Dy();

        // --------- J_E ----------------------------------------------
		double vvx, vvy, vvz, tmp, Ln_ei;

        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;
				//
				// Ni
				tmp = VD.By()(xid,yid).real();

				//
				// grad.Pe
				VD.Ex()(xid,yid).imag() = VD.Bx()(xid,yid).real();
				VD.Ey()(xid,yid).imag() = VD.Bz()(xid,yid).real();
				
                vvx = Vx(ix,iy)        /tmp;
                vvy = Vyz(ix,iy).real()/tmp;
                vvz =-Vyz(ix,iy).imag()/tmp;

				//
				// J/Ni, -grad.Pe/Ni
                VD.Ex()(xid,yid)        /= tmp;
                VD.Ey()(xid,yid)        /= tmp;
                VD.Ez()(xid,yid).real() /= tmp;

				//
				// -Vi(<v>/n) x B + (-grad.Pe/Ni)
                //EF.Ex()(xid,yid).real()  = (vvz*EF.By()(xid,yid).real()-vvy*EF.Bz()(xid,yid).real());
                //EF.Ey()(xid,yid).real()  = (vvx*EF.Bz()(xid,yid).real()-vvz*EF.Bx()(xid,yid).real());
                //EF.Ez()(xid,yid).real()  = (vvy*EF.Bx()(xid,yid).real()-vvx*EF.By()(xid,yid).real());

                EF.Ex()(xid,yid).real()  = (vvz*EF.By()(xid,yid).real()-vvy*EF.Bz()(xid,yid).real()) + VD.Ex()(xid,yid).imag();
                EF.Ey()(xid,yid).real()  = (vvx*EF.Bz()(xid,yid).real()-vvz*EF.Bx()(xid,yid).real()) + VD.Ey()(xid,yid).imag();
                EF.Ez()(xid,yid).real()  = (vvy*EF.Bx()(xid,yid).real()-vvx*EF.By()(xid,yid).real());

				//
				// Dissipation
				//Ln_ei = Eta_ei(tmp, Vsq(ix,iy));
				
                //EF.Ex()(xid,yid).real()  += VD.Ex()(xid,yid).real()*Ln_ei;
                //EF.Ey()(xid,yid).real()  += VD.Ey()(xid,yid).real()*Ln_ei;
                //EF.Ez()(xid,yid).real()  += VD.Ez()(xid,yid).real()*Ln_ei;

            }
        }

        //int moduloX(RANKX()%2); 
		//if(NODESX() >1){
		//	if (moduloX != 0 && Inputdata::IN().list().bndX==2){
		//	}
		//}

			//for(int ix(szx+Nbc-5); ix < szx+Nbc; ++ix){
				//cout<<ix<<": "<<Y.SH(1,0)(2,ix,0)<<endl;
				//cout<<ix<<" Vl: "<<VD.Ex()(ix,3)<<" vs "<<EF.Ex()(ix,3)<<endl;
			//	cout<<ix<<" Vl: "<<VD.Ex()(ix,3)<<" vs "<<Pp(ix,3)<<endl;
			//}
		//
        return;
    }

//--------------------------------------------------------------
// add hall term to E
//--------------------------------------------------------------
    void OHMic::E_Hall(Stat& Yi, Stat& Ye){
		EMF2D& EF(Yi.EMF()), VD(Ye.EMF());
		double vx, vy, vz, tmp;
		int yid, xid;

        for (int iy(0); iy < szy; ++iy){
			yid = iy+Nbc-1;
            for (int ix(0); ix < szx; ++ix){
				xid = ix+Nbc-1;

				// J/Ni
                vx = VD.Ex()(xid,yid).real();
                vy = VD.Ey()(xid,yid).real();
                vz = VD.Ez()(xid,yid).real();

                EF.Ex()(xid,yid).real()  += (vy*EF.Bz()(xid,yid).real()-vz*EF.By()(xid,yid).real());
                EF.Ey()(xid,yid).real()  += (vz*EF.Bx()(xid,yid).real()-vx*EF.Bz()(xid,yid).real());
                EF.Ez()(xid,yid).real()  += (vx*EF.By()(xid,yid).real()-vy*EF.Bx()(xid,yid).real());
            }
        }
        return;
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
// filter routine
//--------------------------------------------------------------
    void OHMic::AUTO_LOW_PC(Stat& Yslope){
       using Inputdata::IN; 
/*
       GSlice_iter< complex<double> > it1(Yslope.SH(1,0).p0(0)), it2(Yslope.SH(1,0).p0(0)),
									  itl(Yslope.SH(1,0).p0(0));
       size_t l0(Yslope.DF().l0()), m0(Yslope.DF().m0());
       size_t maxpeq0(10), tmp_l;
	   //complex<double> A,B;
	   double Ar,Br, Ai, Bi;

	  //--- f(p0) ~p^l
	  //size_t l = 1;{
      for (size_t l = 1; l < l0+1; ++l){
           for (size_t m=0; m<((m0<l)?m0:l+1); ++m){
				for (size_t np=l-1; np <=l;++np){
					if(np >0){
						it2 = Yslope.SH(l,m).p0(l);
						for(itl = Yslope.SH(l,m).p0(np); itl != itl.end(); ++itl)
							*itl = (*it2++) * pow(pr(np)/pr(l+1),l);
					}
				}
				//it2 = Yslope.SH(l,m).p0(1);
				//for(itl = Yslope.SH(l,m).p0(0); itl != itl.end(); ++itl)
				//	*itl = (*it2++) * pow(pr(0)/pr(1),l);
		   }
	  }
	  */
    }
//--------------------------------------------------------------
//
//-------------------------------------------------------------------
     double Eta_ei(double ne, double Te) {
        static double lnei;

		double tmp = 9.627e-39/sqrt(ne*Te*Te*Te)/(Inputdata::IN().list().Zeta*Inputdata::IN().list().Zeta);

        if (ne > 0.0000001) {
            Te /= (3.0*ne);
            Te *= 9.395207885e+08; // assume Te = Ti
            ne *= (Inputdata::IN().list().density_np);

            Te = log(Te); 
            ne = log(ne);

            lnei = 24.0 - 0.5*ne + Te;
			
			//--- e-i
            if (lnei > 2.0){
				lnei *= tmp;

				return lnei;
			}
        }
 
        return  2.0 * tmp;
    }
//-------------------------------------------------------------------
//**************************************************************
