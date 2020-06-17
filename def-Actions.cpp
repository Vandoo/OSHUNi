///////////////////////////////////////////////////////////
//   
//   This cpp file contains the definition for the 
//   collection of methods used for the right hand
//   side of the Runge-Kutta.
///////////////////////////////////////////////////////////
//
//   class Actions::
//
//      The Actions class essentially defines the modules
//      that will be used in the calculation of the right
//      hand side of the Runge-Kutta in the correct order.
//      1. Electric Field
//      2. Current
//      3. Advection
//      4. Magnetic Field
//      5. Faraday's law
//      6. Ampere's law
//      7. A low-pass filter, which is currently deactivated
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
    #include "decl-actions.h"


//**************************************************************
//**************************************************************
//   Definition for the actions
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    //Actions::Actions(Stat& Yslope)
    Actions::Actions(Stat& Yslope, Parallel_Environment& PE)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),  
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0),
         //--- e-i
         E_Act(Yslope, PE),
		 CX(PE),
         //------
         Ef(Yslope),
         Bf(Yslope),
         Jf(Yslope),
         Af(Yslope),
         ME(Yslope){
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Actions::Adv(Stat& Yin, Stat& Ye){ 
//--------------------------------------------------------------
//  Collection of operations
//--------------------------------------------------------------
        Yh  = 0.0;

        //--- e-i
        // explicit solution
		// Yh for temp Ev, New E will be contaminated in RK2/4, but it is updated here again
        //E_Act.Adv_Ohm(Yin, Ye);
        E_Act.Update_E(Yin, Ye);
        //------

        Yh  = Ef.Exyz(Yin);        // Find how E_old changes f !WHY- now it is E_new 

        // This yields the new E
        //------
        //Yh  = Jf.Jxyz(Yin);        // Find the first half of E_new: dE/dt -= J(f_old)  
        //Yh  = ME.Ampere(Yin);      // Find the second half of E_new: dE/dt += rot(B_old)
        //------
        E_Act.Adv_Ohm(Yin, Ye);

        //------
        // This is for either explicit or implicit methods
        Yh  = Af.Axy(Yin);         // Find how A_old changes f  
//        Yh  = Bf.Bxyz(Yin);        // Find how B_old changes f  
//        Yh  = ME.Faraday(Yin);     // Find dB/dt += -rot(E_old)

        CLEAN_LOW_PC(Yh);
        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    //Implicit_Actions::Implicit_Actions(Stat& Yslope)
    Implicit_Actions::Implicit_Actions(Stat& Yslope, Parallel_Environment& PE)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),  
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0),    //WHY? not just use Y.lmax, etc
         if_implicitES(Inputdata::IN().list().if_implicitES),
         if_implicit1D(Inputdata::IN().list().if_implicit1D),
		 CX(PE),
         Ef(Yslope),
         Bf(Yslope),
         Jf(Yslope),
         Af(Yslope),
         ME(Yslope){
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Implicit_Actions::Adv1(Stat& Yin){ 
//--------------------------------------------------------------
//  Collection of operations
//--------------------------------------------------------------
        Yh  = 0.0;

        // This is for either explicit or implicit methods
        Yh  = Af.Axy(Yin);         // Find how A_old changes f  
        
        if (  !(if_implicitES) ) {     // If not electrostatic
            Yh  = Bf.Bxyz(Yin);        // Find how B_old changes f  
            Yh  = ME.Faraday(Yin);     // Find dB/dt += -rot(E_old)
        }

        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Implicit_Actions::Adv2(Stat& Yin){ 
//--------------------------------------------------------------
//  Use the implicit ilectric field to push f
//--------------------------------------------------------------

        Yh  = 0.0;

            //------                      
            // Yh  = Jf.Jxyz(Yin);        // Find the first half of E_new: dE/dt -= J(f_old)  
            // Yh  = ME.Ampere(Yin);      // Find the second half of E_new: dE/dt += rot(B_old)
            //------

        if (  if_implicit1D ) {          // If 1D 
            Yh  = Ef.Exyz1D(Yin);        // Find how E_new changes f  
        }
        else {
            Yh  = Ef.Exyz(Yin);        // Find how E_new changes f  
        }
       

        CLEAN_LOW_PC(Yh);

        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Implicit_Actions::ImpEx(Stat& Yin){ 
//--------------------------------------------------------------
//  Implicit Electric field in the x-direction 
//--------------------------------------------------------------

        Matrix2D< complex<double> > Exi(Yin.EMF().Ex().matrix());

        Yh  = 0.0;

        if (  if_implicit1D ) {                    // If 1D 
            Yh  = Ef.Implicit_Ex1D(Yin, Exi);      // Find how E_x changes f  
        }
        else {                                     // If not 1D  
            Yh  = Ef.Implicit_Ex(Yin, Exi);        // Find how E_x changes f  
        }

        return Yh;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Stat& Implicit_Actions::ImpEy(Stat& Yin){ 
//--------------------------------------------------------------
//  Implicit Electric field in the x-direction 
//--------------------------------------------------------------
        Matrix2D< complex<double> > Eminus(Yin.EMF().Ey().matrix());
        Matrix2D< complex<double> > Eplus(Yin.EMF().Ey().matrix());

        Yh  = 0.0;

        if (  !(if_implicit1D) ) {                      // If not 1D 
            Yh  = Ef.Implicit_Eyz(Yin, Eminus, Eplus);  // Find how E_y changes f  
        }

        return Yh;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Stat& Implicit_Actions::ImpEz(Stat& Yin){ 
//--------------------------------------------------------------
//  Implicit Electric field in the x-direction 
//--------------------------------------------------------------

        complex<double> ii(0.0,1.0);
        Matrix2D< complex<double> > Eminus(Yin.EMF().Ez().matrix());
                                    Eminus *= (-1.0)*ii; 
        Matrix2D< complex<double> > Eplus(Yin.EMF().Ez().matrix());  
                                    Eplus  *=  ii;       

        Yh  = 0.0;

        if (  !(if_implicit1D) ) {                            // If not 1D 
            Yh  = Ef.Implicit_Eyz(Yin, Eminus, Eplus);        // Find how E_z changes f  
        }

        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    //E_Actions::E_Actions(Stat& Yin)
    E_Actions::E_Actions(Stat& Yin, Parallel_Environment& PE)
       : Yin(Yin),
		 CX(PE),
         OM(Yin){
    }
//--------------------------------------------------------------
    void E_Actions::Adv_Ohm(Stat& Yi, Stat& Ye){ 
	  //e.g., ME.Farady needs clean EMF: ...+=.... But the J_Ampere does not	
      //Ye  = 0.0;

	  //--- filter
	  //OM.AUTO_LOW_PC(Yi);

	  //--- Actually updating Ye.field
	  //OM.J_density(Yi, Ye, 0);	 // density profile
      //OM.J_Ampere(Yi, Ye);      // J = rot(B_old)
	  //--- e_pressure
	  //OM.J_AD(Yi, Ye);

	  //--- start with pure E= -UxB MHD case:
	  //    Yh.B has temp N_i, it is removed in ME.Farady 
      //OM.J_E(Yi, Ye);          // Ve = Zi*(Vi-J/n) --->E

	  //--- Tzou pressure
	  OM.J_Pre(Yi, Ye);

	  // need this?, originally for Yh
      //CLEAN_LOW_PC(Yi);

	  //--- for mirror boundary, Ex, By/z continue
	  //CX.Update_XRedge(Yi);

      return;
    }
//--------------------------------------------------------------
    void E_Actions::Update_E(Stat& Yi, Stat& Ye){ 
      
	  OM.E_Hall(Yi, Ye);      // E += JxB

      return;
    }
//--------------------------------------------------------------
//**************************************************************








//**************************************************************
//--------------------------------------------------------------
// This high-harmonic filter is temporarily deactivated.        //WHY? really? it is in act.adv
//--------------------------------------------------------------
    void CLEAN_LOW_PC(Stat& Yh){ 
//--------------------------------------------------------------
//  Sets the low p_cells for the high harmonics equal to zero.
//  Specifically if  (l+m-1 > #p) then #p = 0.0
//--------------------------------------------------------------
       GSlice_iter< complex<double> > it(Yh.SH(0,0).p0(0));     //WHY? useless?
       size_t l0(Yh.DF().l0()), m0(Yh.DF().m0());              
       size_t maxpeq0(80);                                      //WHY? could only change 80 here?
       for (size_t l = 1; l < l0+1; ++l)
           for (size_t m=0; m<((m0<l)?m0:l)+1; ++m) 
               for (size_t np=0; np < (((l-1)<maxpeq0)?(l-1):maxpeq0);++np)       
                   for (GSlice_iter< complex<double> > it(Yh.SH(l,m).p0(np)); 
                       it!=it.end(); ++it) *it = 0.0;
    }
//--------------------------------------------------------------
//**************************************************************



