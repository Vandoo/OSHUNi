
///////////////////////////////////////////////////////////
//   
//   TEMPORARY INPUT
//   This cpp file contains the definitions for the input 
///////////////////////////////////////////////////////////

//
// 
//   namespace Inputdata::
//    The data in the input deck are collected in the class 
//    "Input". For an object obj of class Input the data can 
//    be accessed as follows: 
//
//--> obj.inp().l0  :: (size_t) the l0 number of harmonics
//       "     .m0  :: (size_t) the m0 number of harmonics
//       "     .x   :: Axis<double> for the spatial "x"
//       "     .pr  :: Axis<double> for the radial momentum
//       "     .CLF :: (double) the "Courant" condition
//
//--> obj.cont().n_out  :: (size_t) the number of output steps
//       "      .tstop  :: (double) the end time of the simulation
//       "      .dt_out :: (double) the time between output 
//
//--> obj.outp().x1 :: Axis<float> for the output x1 axis
//       "      .p1 :: Axis<float> for the output p1 axis
//       "      .p2 :: Axis<float> for the output p2 axis
//       "      .p3 :: Axis<float> for the output p3 axis
//
//   
//    Access to the entire list can be obtain from:
//--> obj.list().* 
//   
//    The input is made available to the rest of the program
//    with the function IN(), so that where "obj" above the
//    one will have to substitute "Inputdata::IN()"
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <math.h>
    #include <float.h>

//  My libraries
    #include "matrices.h"

//  Interface
    #include "decl-input.h"



//**************************************************************
//--------------------------------------------------------------
    Inputdata::input_list::input_list(){
//--------------------------------------------------------------
//  The constructor for the input_structure
//  works as a temporary input deck
//--------------------------------------------------------------

//      Parallel
        NnodesX =4; NnodesY = 1;

//      spherical harmonics (assumed l0 > 3, m0 > 1, l0 >= m0)
        l0 = 32;   m0 = 3; 
		l0e= 1;    m0e= 1;
		p0e  = 1;
		nump = 300; numx_glob = 256; numy_glob = 2;

        xmin = -5.12e02;  //  -400um
        xmax =  5.12e02;  //   400um
        ymin = -1.00e00;  //  -400um
        ymax =  1.00e00;  //   400um
        pmax = 3.26247e-4; pmin = pmax /(2.0*nump);

//      the ratio CLF/Dp 
        clf_dp = 2.0e02;
 
//      Algorithms
        if_implicitES = false;
        if_implicit1D = false;

//      Collisions for f00
        small_dt = 5.0e-02;
        smaller_dt = 0.005;
        NB_algorithms = 6;
        if_tridiagonal = true;
        implicit_E     = false;

//      loop control 
        RKLevel = 2; 
        //n_outsteps = 1;
        //t_stop =  50;
        n_outsteps = 2;
        t_stop =  1.0e-2;

//      Restat files
        restart_time = 0; // start reading at the "restart_time"-th restart file
        n_restarts = 0;   // generate "n_restarts" restart files

//      boundary type (default periodic) 
//      0:periodic, 1:mirror, 2: continuity/smooth 
        bndX = 2; bndY = 0;


//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
        sigma_p = 0.22; center_p = 0.5;
//        sigma_x = 0.24; center_x = xmax/2.0;
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//      Numerics output 
        o_EHist = false; 

//      Output options
        o_Ex =   true; o_Ey  = false; o_Ez = false;
        o_Bx =   false; o_By  = false; o_Bz = false;
        o_Jx =   true; o_Jy  = false; o_Jz = false;
		//N
        o_x1x2 = true; o_pth = false;

//      Output stress energy tensor
        o_G  = false;
        o_Px = false;  o_PxPx = false;
        o_Py = false;  o_PxPy = false; o_PyPy = false;
        o_Pz = false;  o_PxPz = false; o_PyPz = false; o_PzPz = false;

//      Nonrelativistic output 
        o_Vx = true; o_VxVx = false; o_VxVz = false;
        o_Vy = false; o_VxVy = false; o_VyVy = false;
        o_Vz = false; o_VyVz = false; o_VzVz = false;
        o_Vsq = true; o_Qx  = false; o_Qy   = false;

        o_Temperature = false;
        o_Pressure = true;
        o_ND = false;
        o_Nu = false;

//      Output 
        o_p1x1 = true;
        o_p1x1x2 = false; 
        o_p1p2p3 = false;
		o_pmulti = 3;
		o_fsp  = true;

//      number of output cells for cartesian data
        nump1 = 128; nump2 = 16; nump3 = 16; 
        numpx = 1000; 

//      Only generate data
        only_output = false;

//      Executable statements
        if (((NnodesX%2))!=0 && (NnodesX !=1)) 
            std::cout<<"The number of nodes "<<NnodesX<<" is not even\n";
        numx      = numx_glob / NnodesX;
        numx     += 2*RKLevel;        
        numx_glob = NnodesX * (numx-2*RKLevel);

        if (((NnodesY%2)!=0) && (NnodesY !=1)) 
            std::cout<<"The number of nodes "<<NnodesY<<" is not even\n";
        numy      = numy_glob / NnodesY;
        numy     += 2*RKLevel;        
        numy_glob = NnodesY * (numy-2*RKLevel);

//      Electron-ion collision parameters
        lnLambda = 10.3; Zeta = 1.0;

//      Density
        density_np = 1.0e+13;// in 1/cm^3

//      Density
        denormalize_fields = false;
        denormalize_length = false;

//      Phenomenological laser parameters
        polarization_direction = 0; // x = 1, y = 2, z = 3, 0 = nothing;
        inverse_bremsstrahlung = 0; // 
        linear_Ivst = false;
        polynomial_Ivst = false;
        spot_w0 = 99999725.41; cntr_r0 = 0.0; cntr_x0 = 0; // -100um 
        I_0 = 1.0e+16; lambda_0 = 0.351; //intensity in W/cm^2, laser wavelength in um 
        rise_time = 0.000000001;//1.7835246; //  
        flat_time = 358927.4021; // 300 ps 
        fall_time = 178352.46; //  100 ps
        ab_eta = 1.0; ab_depth = 5640.0; //1414um 

//---   e-i      
        ex_electric   = false;
        ex_magnetic   = false;
		initial_x_vel = false;
		initial_y_vel = false;
		initial_x_sh  = true;

		Pe_ratio      = 42.850074789; //p_Te/p_Ti
		T_gamma		  = 0.333333333333; //T=n^(2*gamma), gamma = (classic gamma -1)/2
		//constant
		mass_ratio    = 1836.125321096889;
//-----------------------------------------------------------------------
//      Density profile

        // Polynomial
        setup_poly_x = false;
        nmin_x = 0.01;  nconp = 0.01;
        start_x0 = -10000.0, rise_x = 100000.0; flat_x = 120000.0; fall_x = 130000.0;

        // Exponential
        setup_exp_x = false;
        exp_nmin = 1.0e+21;// in 1/cm^3
        exp_xmin = -6.48;// -163um. This is the distance behind 0 which is
                                // where the critical density is, this is determined by the wavelength
                                // The maximum density is density_np

        // Gaussian
        setup_gauss_dens_x = false;
        //xmin = -200.0;   xmax = 512.5;
        ampl_x = 1.0; sigma_x = 8000.0; center_x = 0.0;

        // Sinusoidal density
        setup_sine_dens_x = false;

        // X Piecewise Linear
        setup_pcwise_lin_x = true;
        xloc.push_back( -5.12e02 );   xdens.push_back(   1.0  ); 
        xloc.push_back( -1.00e-01 );   xdens.push_back(   1.0  );
        xloc.push_back(  0.0      );   xdens.push_back(   1.0  );
        xloc.push_back(  1.00e-03 );   xdens.push_back(   1.0  ); 
        xloc.push_back(  2.00e-01 );   xdens.push_back(   1.0  ); 
        xloc.push_back(  5.12e02 );   xdens.push_back(   1.0  );
        //------> Need to check that xloc(0) < xloc(1) < ... etc

        xmin = xloc.at(0); 
        xmax = xloc.at(xloc.size()-1);
        if (xloc.size() != xdens.size()){ 
            std::cout<<"The density profile in 'x' was not declared properly\n";
            exit(1);
        }
        
        // Smoothing level
        smooth_x = 5;

        // Y Gaussian
        setup_gauss_dens_y = false;
        ampl_y = 1.0; sigma_y = 8000.0; center_y = 0.0;
        //ymin = -570.0; ymax = 570.0;

        // Y Piecewise Linear
        setup_pcwise_lin_y = true;
        yloc.push_back( -1.00e00);   ydens.push_back(   1.0  ); 
        yloc.push_back( -8.00e-04);   ydens.push_back(   1.0  );
        yloc.push_back( -4.00e-04);   ydens.push_back(   1.0  );
        yloc.push_back(  6.00e-04);   ydens.push_back(   1.0  ); 
        yloc.push_back(  8.00e-04);   ydens.push_back(   1.0  ); 
        yloc.push_back(  1.00e00);   ydens.push_back(   1.0  );
        //------> Need to check that yloc(0) < yloc(1) < ... etc

        ymin = yloc.at(0); 
        ymax = yloc.at(yloc.size()-1);
        if (yloc.size() != ydens.size()){ 
            std::cout<<"The density profile in 'y' was not declared properly\n";
            exit(1);
        }
 //-----------------------------------------------------------------------

 //-----------------------------------------------------------------------
 //     Temperature profile

        // Minimum / ambient temperature
        pt_amb = 3.262472e-05;

        // Gaussian
        setup_gauss_temp_x = false;
        sigma_tx = 16000.0; center_tx = 0.0; ptx = 4.51e-2;

        // Sinusoidal Temperature
        setup_sine_temp_x = false;
        pt_sinmax = 1.0;

        // X
        setup_pcwise_tmp_x = true;
        xtloc.push_back(  -5.12e02 );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(  -1.32e-03 );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(  -0.48e-04 );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(   0.0      );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(   1.32e-04 );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(   1.24e-03 );   xtemp.push_back(3.262472e-05);
        xtloc.push_back(   5.12e02 );   xtemp.push_back(3.262472e-05);


        if (xtloc.size() != xtemp.size()){ 
            std::cout<<"The temperature profile in 'x' was not declared properly\n";
            exit(1);
        }

        // Gaussian
        setup_gauss_temp_y = false;
        sigma_ty = 16000.0; center_ty = 0.0; pty = 1.0;
        
        // Sinusoidal Temperature
        setup_sine_temp_y = false;
        pt_sinmay = 3.4256e-5;

        // Y
        setup_pcwise_tmp_y = true;
        ytloc.push_back(  -1.00e00 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(  -1.72e-04 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(  -2.99e-05 );   ytemp.push_back(   1.00  );
        ytloc.push_back(   2.99e-05 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(   1.72e-04 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(   1.00e00 );   ytemp.push_back(   1.00  );

        if (ytloc.size() != ytemp.size()){ 
            std::cout<<"The temperature profile in 'y' was not declared properly\n";
            exit(1);
        }
 //-----------------------------------------------------------------------

        return;
    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------

//    Definition of the constructor for the input data struct
      Inputdata::decl_input:: decl_input(input_list in) :
          l0(in.l0), m0(in.m0), 
          pr(in.nump, in.pmin, in.pmax),
          x( in.numx, in.xmin, in.xmax), 
          xglob(in.numx_glob, in.xmin, in.xmax), 
          y( in.numy, in.ymin, in.ymax), 
          yglob(in.numy_glob, in.ymin, in.ymax), 
          CLF(/*pr.dx()**/in.clf_dp) { 

          x(0) = xglob(0)-in.RKLevel*xglob.dx();
          for (int i = 1; i< x.dim(); ++i) x(i) = x(i-1) + xglob.dx();
          y(0) = yglob(0)-in.RKLevel*yglob.dx();
          for (int i = 1; i< y.dim(); ++i) y(i) = y(i-1) + yglob.dx();
          return; 
       }
           

//--------------------------------------------------------------
//    Definition of the constructor for the loop control data struct
      Inputdata::decl_control::decl_control(input_list in) :
          n_out(in.n_outsteps), tstop(in.t_stop)
          { dt_out = tstop/double(n_out); return; }
//--------------------------------------------------------------

//--------------------------------------------------------------

//    Definition of the constructor for the output data struct
      Inputdata::decl_output::decl_output(input_list in) :
 //         x1(in.numx1,static_cast<float>(in.xmin), static_cast<float>(in.xmax)), 
          px(in.numpx,static_cast<float>(in.pmax)), 
          p1(in.nump1,static_cast<float>(in.pmax)), 
          p2(in.nump2,static_cast<float>(in.pmax)), 
          p3(in.nump3,static_cast<float>(in.pmax)) 
          { return; }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the Input Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    Inputdata::Input:: Input(){
        ilist = new input_list();
        in = new decl_input(*ilist); 
        control = new decl_control(*ilist); 
        out = new decl_output(*ilist);
    }
//  Destructor
    Inputdata::Input:: ~Input(){
        delete ilist;
        delete in;
        delete control;
        delete out;
    }

//--------------------------------------------------------------
//  Access 
//--------------------------------------------------------------
    Inputdata:: input_list& Inputdata::Input:: list() const  {return (*ilist);}
    Inputdata:: decl_input& Inputdata::Input:: inp() const {return (*in);}
    Inputdata:: decl_control& Inputdata::Input:: cont() const {return (*control);}
    Inputdata:: decl_output& Inputdata::Input:: outp() const {return (*out);}

//--------------------------------------------------------------
//**************************************************************


//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
    Inputdata::Input& Inputdata::IN() {
        static Inputdata::Input in;
        return in;
    }

//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
