///////////////////////////////////////////////////////////
//
//   class Actions::
//
//      The Actions class essentially defines the modules
//      that will be used in the calculation of the right
//      hand side of the Runge-Kutta.
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

    #ifndef DECLERATION_ACTIONS_H
    #define DECLERATION_ACTIONS_H

//**************************************************************
//**************************************************************
//   Definition for the Output namespace
//**************************************************************
//**************************************************************

#include "decl-export.h"
#include "decl-parallel.h"


//**************************************************************
        void CLEAN_LOW_PC(Stat& Yh);
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class E_Actions {
        private:
            Stat& Yin;
            OHMic OM;
			Parallel_Environment& CX;
        public:
//          Constructors/Destructors
            //E_Actions(Stat& Yin); 
            E_Actions(Stat& Yin, Parallel_Environment& PE); 
         
//          Explicit Advance
            void Adv_Ohm(Stat& Yi, Stat& Ye); 
            void Update_E(Stat& Yi, Stat& Ye); 
		};
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Actions {
//--------------------------------------------------------------
//      Decleration of the Actions
//--------------------------------------------------------------
        private:
//          Define the Actions
            Stat& Yh;
            ElectricField Ef;
            MagneticField Bf;
            Current Jf;
            Spatial_Advection Af;
            MaxwellEq ME;
            //--- e-i
            E_Actions E_Act;
			Parallel_Environment& CX;

            size_t l0, m0;

        public:
//          Constructors/Destructors
            //Actions(Stat& Yslope);
            Actions(Stat& Yslope, Parallel_Environment& PE);
         
//          Explicit Advance
            Stat& Adv(Stat& Yin, Stat& Ye);

        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        class Implicit_Actions {
//--------------------------------------------------------------
//      Decleration of the Actions
//--------------------------------------------------------------
        private:
//          Define the Actions
            Stat& Yh;
            ElectricField Ef;
            MagneticField Bf;
            Current Jf;
            Spatial_Advection Af;
            MaxwellEq ME;

            bool implicit_E;
            bool if_implicitES;
            bool if_implicit1D;
            size_t l0, m0;

        public:
//          Constructors/Destructors
            //Implicit_Actions(Stat& Yslope); 
            Implicit_Actions(Stat& Yslope, Parallel_Environment& PE); 
         
//          Explicit Advance
            Stat& Adv1(Stat& Yin);    
            Stat& Adv2(Stat& Yin);    

//          Implicit Advance
            Stat& ImpEx(Stat& Yin);    
            Stat& ImpEy(Stat& Yin);    
            Stat& ImpEz(Stat& Yin);    
			// --- e-i
			Parallel_Environment& CX;
        };
//--------------------------------------------------------------
//**************************************************************

    #endif
