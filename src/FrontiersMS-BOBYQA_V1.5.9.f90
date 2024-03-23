! 31-Dec-2017 V1.5.9
!  -  Increased sharpness of step up for delG for Bac and SRB slightly
!  -  Uses new function fZero to avoid state variables going less than or equal to zero.
! 28-Dec-2017 V1.5.8
!  -  Decomposition reactions for Bac and SRB only occur if delG < 0. 
! 22-Dec-2017 V1.5.7
!  -  Added kIntg_x to the light output file (envUnit)
! 20-Dec-2017 V1.5.6
!  -  Added output for sigmaDot for each reaction, and also added sumSigmaDot for reactions, water
!     and particles to output.
! 17-Dec-2017 V1.5.5
!  -  Noticed that weightS was not being applied to sumSigma. This was corrected.
! 16-Dec-2017 V1.5.4
!  -  findMaxEntropy did not end with a call to entropy, so sumSigmaTotInf was not be returned
!     with the correct value, and explains why sumSigmaTotInf differs from sumSigmaTotFin with the
!     when the two intervals are the same.
! 11-Dec-2017 V1.5.3
!  -  In integrateState, sumSigmaDotPrev was set to zero at the start of an integration, but that's
!     not correct, as sumSigmaDotPrev can be determined from the initial conditions at t0. 
!     This was corrected, but it should not make much difference.
!  9-Dec-2017 V1.5.2
!  -  Changed initialization to read uppercase omega instead of lower case omega.  This way,
!     the provious solution can be used to inialize a continuation run.
!  5-Dec-2017 V1.5.1
!  -  Added errCnt to screen output
!  -  Added output of PAR in uE/m2/s to be more comparible to typical units (6-Dec)
!  3-Dec-2017 V1.5
!  -  Added KCOL to the parameter input list to allow changes to that KCOL
!  -  Also, added better stopping execution if any process has an integration failure during
!     solution saving.
!  -  Added reporting myRank in integration failure
!  2-Dec-2017 V1.4.2
!  -  Initializing the PDE is still a problem, so also try introducing spinup of q and qL
! 26-Nov-2017 V1.4.1
!  -  Added MPI_Barrier to findMaxEntropy routine.  I don't think it's necessary but it does not hurt
!     in this case.
! 15-Nov-2017 V1.4
!  -  Replaced nTpts with deltaT, where deltaT is the time spacing, in days, between PDE integration
!     time points. This way, integration over deltInterval has the same spacing as over deltInfinity,
!     or at least very close (note, deltaT may not be the exact value, as it will be adjusted so that
!     the PDE integrator will integrate exactly to tf.
! 30-Oct-2017 V1.3
!  -  This is updated from V.17 of pCMA and 1.2 of vtDirect
!  -  4-Nov-2017: made some minor improvements in displaying time outputs
! 29-Sep-2017 V1.2
!  -  This version for FrontiersMS that uses BOBYQA was derived from V1.5 of pCMA version, as V1.5
!     containts the latest biogeochemistry and other programming modifications.  This version als
!     replaces OpenMP with MPI in the main loop that optimizes over depth. OpenMP is used, as in V1.5
!     to speed up BALCOLI.  
   module realPrec
      ! this sets the precision of reals
      implicit none
      save
      integer, parameter:: sp = kind(1.0E0)
      integer, parameter:: dp = kind(1.0D0)
   end module realPrec
   
   module interfaces
      interface
         subroutine omega_free (n, w, omg)
            use realPrec
            implicit none
            integer, intent(in)  :: n              ! number of sub reactions for a biological structure
            real(dp), intent(in) :: w(*)  ! lower case omega, control variables
            real(dp), intent(out):: omg(*)         ! upper case omega, these sum to 1, which is where the nth omega comes from
         end subroutine omega_free
      end interface
      interface
         pure subroutine insertion_sort(n, a)
            use realPrec
            implicit none
            integer, intent(in)                    :: n
            real(dp), intent(in out), dimension(:) :: a
         end subroutine insertion_sort
      end interface   
   end module interfaces
   
   module typeVars
      ! This defines the derived type to pass variables and parameters to
      ! the subroutines that calculate reaction rates, etc.
      ! This allows free energies of formation, etc to be calculated once, then
      ! passed to routine to allow for thread safe operation.
      ! If a new reaction needs a new parameter, the rxnVariables derived type
      ! can just be added to.  
      !
      ! This routine is problem specific, so must be changed as reactions are changed, etc
      use realPrec
      implicit none
      type rxnVariables
         real(dp) T_K         ! Local temperature (K) at time and location x
         real(dp) csAtx       ! local cross sectional area (m2)
         real(dp) pH          ! pH at depth x
         real(dp) is          ! local ionic strenth (M) at (t,x)
         real(dp) depth       ! water column depth (m)
         real(dp) absZero     ! what is considered 0., but not equal to that to avoid division by 0.
         real(dp) dGf0_bioS   ! free energies of formation of species at current temp (K), is (M), and pH. (kJ/mol)
         real(dp) dGf0_C_D    ! Detritus C, N and P.
         real(dp) dGf0_N_D
         real(dp) dGf0_P_D
         real(dp) dGf0_ch2o   ! glucose, single carbon
         real(dp) dGf0_h2aq   ! this are calculated using thermoData in routine rxnProperties 
         real(dp) dGf0_o2aq
         real(dp) dGf0_h2so4
         real(dp) dGf0_h2saq
         real(dp) dGf0_dic
         real(dp) dGf0_hno3
         real(dp) dGf0_nh3
         real(dp) dGf0_h3po4
         real(dp) dGf0_h2o
         real(dp) dGr_Ggamma  ! Energy available in green light, accounting for thermo efficiency (J/mmol)
         real(dp) co2         ! free co2  concentration at time t and depth x (mmol/m3)
         real(dp) hco3        ! free hco3 concentration at time t and depth x (mmol/m3)
         real(dp) RkJ         ! gas constant in J/K/mmol or kJ/K/mol
         real(dp) nuStar      ! max growth rate (mmol/m3/d)
         real(dp) nuDet       ! max growth rate on detritus (mmol/m3/d)
         real(dp) kappa       ! substrate affinity (mmol/m3)
         real(dp) alp         ! Use for CHONP of biological structures (all the same for now, and single C basis)
         real(dp) bet
         real(dp) gam
         real(dp) del
         real(dp) k_w         ! light absorption by water (1/m)
         real(dp) k_p         ! light absoprtion coef by algae and other particles (m^2/mmol-S)
         real(dp) cell_F      ! The true concentration of C_P is much higher that the state solution for because it is only inside the cell.  Adjust for kinetics
         real(dp) Itx         ! light intensity at the at time t and depth x (mmol photons /m^2 /d (Not micromoles))         
      end type rxnVariables
      
      type concentrations
         ! this is mostly used to pass concentration state variables for biological structure calculations
         ! it includes the biological structure concentrations too.
         ! all in mmol/m3 except sal, which is in PSU
         real(dp) sal   ! 1)
         real(dp) o2    ! 2)
         real(dp) dic   ! 3)
         real(dp) nh3   ! not used as a state var yet
         real(dp) po4   ! 4)
         real(dp) so4   ! 5)
         real(dp) h2s   ! 6)
         real(dp) C_Phy ! 7) Carbon storage in phototrophs
         real(dp) C_GSB ! 8)  "      "      "  Green sulfur bacteria
         real(dp) C_L   ! 9) Dissolved labile carbon
         real(dp) C_D   ! 10) Detrital C
         real(dp) N_D   !   "      N BUT not used yet as a state var yet
         real(dp) P_D   ! 11)  "      P
         real(dp) PhyS  ! 12) phytoplankton biological structure
         real(dp) GSBS  ! 13) biological structure
         real(dp) GzS   ! 14) biological structure
         real(dp) AGzS  ! 15) biological structure
         real(dp) BacS  ! 16) biological structure
         real(dp) SRBS  ! 17) biological structure
         real(dp) PhS   ! 18) biological structure
         real(dp) SOxS  ! 19) biological structure
      end type concentrations
      
      type PhyCntl
         ! Conctrol variables associated with Phy biological structure
         real(dp) eps
         real(dp) omg(2) ! note, these are the capital omegas, calculated from w1 control variable
      end type PhyCntl            
      type PhyRxns
         ! These are output variables for the Phy reactions
         real(dp) r(2)   ! reaction rates mmol/m3/d
         real(dp) dG(2)  ! reaction free energies J/mmol or kJ/mol
         real(dp) ne(2)  ! Number of electrons transfered per reaction
         real(dp) sig(2) ! internal entropy production for each reaction (J/K/d/m)
         real(dp) n(2)   ! catabolic rxn driver to anabolic rxn so delrG = 0 when eps = 1 (dimensionless)
         real(dp) aA2, bA2 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(2), F_K(2) ! Thermodynamic and kinetic reaction drivers
      end type PhyRxns
      
      type GSBCntl
         ! Control variables associated with GSB (green sulfur bacteria)
         real(dp) eps
         real(dp) omg(2) ! note, these are the capital omegas, calculated from w1 control variable
      end type GSBCntl               
      type GSBRxns
         ! Outputs associated with GSB reactions
         real(dp) r(2)   ! reaction rates mmol/m3/d
         real(dp) dG(2)  ! reaction free energies J/mmol or kJ/mol
         real(dp) ne(2)  ! Number of electrons transfered per reaction
         real(dp) sig(2) ! internal entropy production for each reaction (J/K/d/m)
         real(dp) n(2)   ! catabolic rxn driver to anabolic rxn so delrG = 0 when eps = 1 (dimensionless)
         real(dp) aA2, bA2 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(2), F_K(2) ! Thermodynamic and kinetic reaction drivers
      end type GSBRxns    

      type GzCntl
         ! Control variables associated with Gz (grazers)
         real(dp) eps
         real(dp) omg(8) ! note, these are the capital omegas, calculated from w1 control variable
      end type GzCntl               
      type GzRxns
         ! Outputs associated with Gz reactions (a total of 7)
         real(dp) r(8)   ! reaction rates mmol/m3/d
         real(dp) dG(8)  ! reaction free energies J/mmol or kJ/mol
         real(dp) nei    ! Number of electrons transfered per reaction (all the same)
         real(dp) sig(8) ! internal entropy production for each reaction (J/K/d/m)
         real(dp) aCi, bCi ! coefficicents to balance Catabolic reaction, but they are all the same.
         real(dp) F_T(8), F_K(8) ! Thermodynamic and kinetic reaction drivers
      end type GzRxns    
      
      type AGzCntl
         ! Control variables associated with AGz (anerobic grazers)
         real(dp) eps
         real(dp) omg(8) ! note, these are the capital omegas, calculated from w1 control variable
      end type AGzCntl               
      type AGzRxns
         ! Outputs associated with AGz reactions (a total of 7)
         real(dp) r(8)   ! reaction rates mmol/m3/d
         real(dp) dG(8)  ! reaction free energies J/mmol or kJ/mol
         real(dp) nei    ! Number of electrons transfered per reaction (all the same)
         real(dp) sig(8) ! internal entropy production for each reaction (J/K/d/m)
         real(dp) aCi, bCi ! coefficicents to balance Catabolic reaction, but they are all the same.
         real(dp) F_T(8), F_K(8) ! Thermodynamic and kinetic reaction drivers
      end type AGzRxns  
      
      type BacCntl
         ! Control variables associated with Bac (heterotrophic bacteria)
         real(dp) eps
         real(dp) omg(3) ! note, these are the capital omegas, calculated from w1 control variable
      end type BacCntl               
      type BacRxns
         ! Outputs associated with Bac reactions
         real(dp) r(3)     ! reaction rates mmol/m3/d
         real(dp) dG(3)    ! reaction free energies J/mmol or kJ/mol
         real(dp) ne       ! Number of electrons transfered per reaction
         real(dp) sig(3)   ! internal entropy production for each reaction (J/K/d/m)
         real(dp) aA1, bA1 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(3), F_K(3) ! Thermodynamic and kinetic reaction drivers
      end type BacRxns        
  
      type SRBCntl
         ! Control variables associated with SRB (sulfate reducing bacteria)
         real(dp) eps
         real(dp) omg(3) ! note, these are the capital omegas, calculated from w1 control variable
      end type SRBCntl               
      type SRBRxns
         ! Outputs associated with SRB reactions
         real(dp) r(3)     ! reaction rates mmol/m3/d
         real(dp) dG(3)    ! reaction free energies J/mmol or kJ/mol
         real(dp) ne       ! Number of electrons transfered per reaction
         real(dp) sig(3)   ! internal entropy production for each reaction (J/K/d/m)
         real(dp) n1       ! catabolic rxn driver to anabolic rxn so delrG = 0 when eps = 1 (dimensionless)
         real(dp) aA1, bA1 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(3), F_K(3) ! Thermodynamic and kinetic reaction drivers
      end type SRBRxns   
      
      type PhCntl
         ! Control variables associated with Ph (photoheterotrophs)
         real(dp) eps
      end type PhCntl               
      type PhRxns
         ! Outputs associated with SRB reactions
         real(dp) r(1)     ! reaction rates mmol/m3/d
         real(dp) dG(1)    ! reaction free energies J/mmol or kJ/mol
         real(dp) ne       ! Number of electrons transfered per reaction
         real(dp) sig(1)   ! internal entropy production for each reaction (J/K/d/m)
         real(dp) aA1, bA1 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(1), F_K(1) ! Thermodynamic and kinetic reaction drivers
      end type PhRxns      
      
      type SOxCntl
         ! Control variables associated with SOx (sulfur bacteria)
         real(dp) eps
      end type SOxCntl               
      type SOxRxns
         ! Sulfide oxidizers (sulfur bacteria)
         real(dp) r(1)     ! reaction rates mmol/m3/d
         real(dp) dG(1)    ! reaction free energies J/mmol or kJ/mol
         real(dp) ne       ! Number of electrons transfered per reaction
         real(dp) sig(1)   ! internal entropy production for each reaction (J/K/d/m)
         real(dp) n1       ! catabolic rxn driver to anabolic rxn so delrG = 0 when eps = 1 (dimensionless)
         real(dp) aA1, bA1 ! coefficicents to balance Anabolic and/or Catabolic reactions
         real(dp) F_T(1), F_K(1) ! Thermodynamic and kinetic reaction drivers
      end type SOxRxns      
      
      type rxnNet
         ! this allows all the calculated values from each the biological structure to be 
         ! passed back to the calling rouine needed for further calculations.
         ! This type now also holds the values for the control variables at a given time and depth as well.
         type(PhyCntl) PhyC ! Phy Control variables
         type(PhyRxns) PhyR ! Phy reaction rates, stoichiometry and thermodynamics
         type(GSBCntl) GSBC ! GSB control variables
         type(GSBRxns) GSBR ! GSB reaction rates, stoichiometry and thermodynamics
         type(GzCntl ) GzC  ! Gz control variables
         type(GzRxns ) GzR  ! Gz reaction rates, stoichiometry and thermodynamics
         type(AGzCntl) AGzC ! AGz control variables
         type(AGzRxns) AGzR ! AGz reaction rates, stoichiometry and thermodynamics
         type(BacCntl) BacC ! Bac control variables
         type(BacRxns) BacR ! Bac reaction rates, stoichiometry and thermodynamics
         type(SRBCntl) SRBC ! SRB control variables
         type(SRBRxns) SRBR ! SRB reaction rates, stoichiometry and thermodynamics
         type(PhCntl ) PhC  ! Ph control variables
         type(PhRxns ) PhR  ! Ph reaction rates, stoichiometry and thermodynamics
         type(SOxCntl) SOxC ! SOx control variables
         type(SOxRxns) SOxR ! SOx reaction rates, stoichiometry and thermodynamics
      end type rxnNet     
      ! It turns out the allocatable and automatic arrays is subroutines can take a crap load of time in allocation/dellalocation
      ! Since the parameters below are problems specific (they can't be changed without modifying the code, it much faster to
      ! set them as parameters and use them when declaring arrays.  This will result in more memory use, but that's OK.
      integer, parameter:: maxSubRxn = 8 ! This is the maximum number of sub reactions in any biological structure
      integer, parameter:: maxPDE = 19 ! the number of PDE state variables (nBioS+nConc).  
      ! This is the number of subreactions for each biological structure in the order they are stored in bioS
      ! In version 1.5.2 and greather, this is used to read in cntl variables that were stored with OMEGA and not w.
      integer, parameter:: nSubRxns(8) = [2, 2, 8, 8, 3, 3, 1, 1]
   end module typeVars    
   
   module parameters
      use realPrec
      use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown       
      implicit none
      save
      integer nconc ! number of concentration variables
      integer ncntl(3) ! number of time control points, number of space control points, total number of control variables including ohmegas
      integer nbioS ! number of biological structures
      ! Parameters associated with problem dimensions
      namelist /dimens/ nconc, ncntl, nbioS
      
      ! Parameters associated with BACOLI95
      integer nXpts ! number of evenly space x points to save at time t.  Note, spacial domain should be normalized between 0 and 1.
      integer nXpts_max ! maximum number of interior points BACOLI is allowed to use
      !integer nTpts ! Number of time points to save for an interval in the PDE soln. nTpts = 1 means just save tf     
      real(dp) deltaT ! time spacing (days) between time outputs in PDE integration
      real(dp) atol1_bacol(1), rtol1_bacol(1) ! scalers for absolute and relative toldrances for BACOLI95
      integer kcol_bacol ! parameter value for kcol. This is set in bacoli95_init
      namelist /BACOLIparams/ nXpts, nXpts_max, deltaT, atol1_bacol, rtol1_bacol, kcol_bacol
      real(dp), allocatable:: Xpts(:) ! Xpts to get solution that will be stored in Cout, along with first spatial derivatives
      logical saveSoln_BACOLI ! whether BACOLI saves solution or not.
      
      type(bacoli95_sol):: ADRsol ! holds PDE solution from bacoli95
      logical solGood ! is true if a good solution exist in ADRsol
      real(dp), allocatable:: cntl(:,:,:) ! number of control variables by time by space, BUT, this should not be placed in a module, for later openMP development.
      real(dp), allocatable:: concIni(:,:) ! Holds the initial conc coniditon over an interval (nconc, nXpts)
      real(dp), allocatable:: bioSini(:,:) ! Holds initial biological structures (nbioS, nXpts)
      real(dp), allocatable:: cntlIni(:,:) ! holds the initial values of control vector over space at t0 also in IC file
      real(dp), allocatable:: concBC(:,:) ! Concentration of conc in feed (nconc dim)
      real(dp), allocatable:: bioSBC(:,:) ! Concentration of bioS in feed (nbioS dim)
      real(dp), allocatable:: tOCnodes(:) ! Vector of time grid points for optimal control (0:ncntl(1))
      real(dp), allocatable:: xOCnodes(:) ! Vector of space grid points for optimal control (ncntl(2))      
      real(dp), allocatable:: lowerBnd(:), upperBnd(:) ! The upper and lower bounds on the control variables for all t and x.
      real(dp), allocatable:: cLat(:)  ! laterial inputs, assumed fixed over time and laterial space.
      real(dp), allocatable:: vSink(:) ! sinking velocity for each state variable read from *.inp (m/d)
      real(dp), allocatable:: kIntg(:) ! This stores the integral of the light abstorption, BUT from the previous saved time point. (a crude approximation)
      integer iNode   ! current depth being optimized

      real(dp) absZero ! Number to use as zero    
      real(dp) depth ! depth of system (m) Note, it is assumed that the over boundary is at 0 m.  For now, use dimensional units (m) 
      real(dp) k_w   ! light absorption by water (1/m)
      real(dp) k_p   ! light absoprtion coef by algae and other particles (m^2/mmol-S)
      real(dp) I0max ! Maximum Light intensity at surface (mmol photons /m^2 /d (Not micromoles/s))
      real(dp) dLat  ! Latitude for calculating solar radition
      real(dp) delrGgamma ! Energy available in green light, accounting for thermo efficiency (J/mmol)
      real(dp) kappa ! half saturation constant
      real(dp) nuStar ! max growth rate
      real(dp) nuDet  ! max growth rate on Detritus
      real(dp) tStart, tFinish ! start and stop time of the interval (d)
      real(dp) deltInfinity ! Length of the infinite interval (d)
      real(dp) deltInterval ! Length of the update interval (d)
      real(dp) wInfinity    ! value of weighting (discounting) function at tInfinity; range: (0,1].  Value of 1 means no weighting
      real(8) WTime_0, WTime_f ! used for timing
      logical useSigmaW  ! if set to true, the total entropy production is weight be weightS (tInt0, t).  Used for optimal control.
      integer maxLocalIter, localIter ! maximum number of times local iteration to find local maximums is allowed, and the actual value used.
      real(dp) maxLocalRes, localRes ! the maximum residual between succesive local max iternation, and the actual local residual
      logical PDEonly ! if set to true, then only run the PDE solution for one infinate interval.
      integer interval ! interval solution is current searching in
      !integer readADRsol ! if >0, then a file basefilename.sol will be read to recover a previously stored BACOLI solution given by the number
      logical readADRsol ! if true then a file basefilename.sol will be read to recover a previously stored BACOLI solution.
      real(dp) pV_o2, pV_co2, pV_h2s ! piston velocities for O2, CO2 and H2S (m/d)
      real(dp) pO2, pCO2, pH2S ! partial pressures of O2 CO2, and H2S in the atmosphere above water (atm)
      real(dp) f_OMremin ! fraction of setteling particulate OM that is remineralized by benthose (unitless)
      real(dp) k_O2remin, k_SO4remin  ! half saturation constant for O2 and SO4 consumption by benthose. Just make it small (uM)
      real(dp) bioS_H, bioS_O, bioS_N, bioS_P ! elemental composition of all biological structures, unit carbon based.
      real(dp) nh3_hold, N_D_hold ! NH3 and detrital N constant conc. used for thermo and kinetic calculations
      real(dp) cell_F ! concentration factor for C_P. Intracellular versus extracellular volume.  In theory this changes as phyS conc does, but keep as constant for now.
      real(dp) delPsi ! this is used in the LaRowe2012 thermodynamic driver function. It is the membrane potential in Volts (not mV)
      real(dp) tStep, tSig ! parameters that allow graceful increase to POM to benthose to help PDE initialization.  tStep is where the flux is 1/2 max, and tSig is the exponent.
      integer errCnt, errCntTot ! Number of errors that occur in an interval in a process, and all the errors.
      integer myRank ! processes rank for MPI.
      integer fcnCalls, fcnCallsTot ! total number of fucntion calls by a process, and in total
      integer ompThreads ! The number of OpenMP threads to use by each MPI process.
      
      ! To remove concetration and space dimensions
      real(dp) Xc ! characterist lenght scale
      real(dp), allocatable:: Cc(:) ! characteristic concentrations and other state variable values

      character (len=80) basefilename
      integer inpUnit, summaryUnit, icUnit, gridUnit, OCgridUnit, solUnit ! units given to input and output streams that are sequential files
      ! Files opened for Direct Access.  Note, these files MUST have the same format.
      integer concUnit
      integer bioSUnit  
      integer cntlUnit  
      integer rxnUnit
      integer dGrUnit
      integer F_Kunit   ! Unit to save F_K reaction driver info
      integer F_Tunit   ! Unit to save F_T reaction driver info
      integer sigmaUnit        
      integer envUnit   ! output for various environmental variables
      integer flxUnit   ! saves the flux of consituents.
      integer BACunit   ! Unit to save BACOLI95 solution info
      integer mpiBOBout ! MPI "unit" for output from BOBYQA
      integer sigmaRijUnit  ! unit that stores entropy production associated with reactions
      integer, parameter:: IJKrec = 3 ! this is the record number where the tecplot IJK info is saved.
      integer  iTec, jTec, kTec ! used for TecPlot
      character(len=2), parameter:: CRLF = char(13)//char(10) ! used by MPI writes for error output
      
      ! The guess as to the numerical precision of calculating total internal entropy production
      ! use 0. if unknown
      real(dp) precSigma
      
      ! Parameters associated with problem
      namelist /params/ absZero, depth, k_w, k_p, I0max, dLat, delrGgamma 
      namelist /params/ kappa, nuStar, nuDet, tStart, tFinish 
      namelist /params/ precSigma, deltInfinity, deltInterval, wInfinity
      namelist /params/ lowerBnd, upperBnd, maxLocalIter, maxLocalRes
      namelist /params/ xOCnodes ! the ncntl(2) nodes in x grid used for the OC grid. They can be spaced however makes sense. tOCnodes are generated by program.
      namelist /params/ concBC, bioSBC 
      namelist /params/ PDEonly
      namelist /params/ Xc, Cc, cLat
      namelist /params/ readADRsol
      namelist /params/ vSink
      namelist /params/ pV_o2, pV_co2, pV_h2s ! piston velocities for O2, CO2 and h2s (m/d)
      namelist /params/ pO2, pCO2, pH2S ! partial pressures of O2, CO2 and hs2 in the atmosphere above water (atm)
      namelist /params/ f_OMremin ! fraction of setteling particulate OM that is remineralized by benthose (unitless)
      namelist /params/ k_O2remin, k_SO4remin  ! half saturation constant for O2 and SO4 consumption by benthose. Just make it small (uM)
      namelist /params/ bioS_H, bioS_O, bioS_N, bioS_P ! elemental composition of phytoplankton, unit carbon based.
      namelist /params/ nh3_hold, N_D_hold ! NH3 and detrital N constant conc. used for thermo and kinetic calculations
      namelist /params/ cell_F ! concentration factor for C_P. Intracellular versus extracellular volume.  In theory this changes as phyS conc does, but keep as constant for now.
      namelist /params/ delPsi ! parameter in LaRowe2012 thermodynamic driver function (volts)
      namelist /params/ tStep, tSig ! parameters that allow graceful increase to POM to benthose to help PDE initialization.  tStep is where the flux is 1/2 max, and tSig is the exponent.
      namelist /params/ ompThreads
      
      ! parameters used by BOBYQA
      integer npt_bobyqa ! Number of interpolation points
      real(dp) rhobeg, rhoend ! initial and final values of a trust region radius
      integer iprint ! controls amount of printing (0, 1, 2 or 3)
      integer maxfun ! maximum number of calls to CALFUN
      namelist /BOBYQAparams/ npt_bobyqa, rhobeg, rhoend, iprint, maxfun  
      
   contains
      subroutine initialize
         ! gets base filename, openFiles for input and output, reads parameters and allocates some space. 
         use mpi
         use mpiErrOut ! this defines mpiFHerr that can be passed to bacoli_errOut_MPI.f
         use typeVars, only: nSubRxns
         implicit none
         integer narg, i, j, ioerr, ioE, mpierr, ncntlOMG, iOMG, iWW
         character (len=2000) string
         character (len=30) tmp
         real(dp), allocatable:: temp(:)
         
         ! First try getting basefilename from commandline
         ! if that fails, then have processes zero get name and broadcast it to other processes 
         narg = command_argument_count () ! see if a filename has been included in command line
         if (narg == 1) then
            call get_command_argument(narg,value=basefilename)
         else 
            if (myRank == 0) then
               write(6,'(a,$)') 'Enter parameter filename base, no extension: '
               read(5,*) basefilename
            end if
            call MPI_BCAST(basefilename,80,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
         end if
         ! open a readonly unit to the inp file. This has to be readonly so that each process can open
         ! it without conflict. This file must be available for all processes.
         open(newunit=inpUnit,file=trim(basefilename)//'.inp',action='read',status='old',iostat=ioerr)   
         if (ioerr /= 0) then 
            write(6,'(a,i5,a)') 'Process ',myRank,' could not open INP file' ! this may or maynot print.
         else
            ! get dimension of problem, and allocate space
            read(inpUnit,nml=dimens)
            ! get BACOLI parameters
            read(inpUnit,nml=BACOLIparams)     
            ! get BOBYQA parameters
            read(inpUnit,nml=BOBYQAparams)   
            if (nbios /= size(nSubRxns)) then
               write(*,'(a)') 'nbios not equal so the size of nSubRxns'
               ioerr = 1
            end if
            ncntlOMG = nbioS + sum(nSubRxns, mask=nSubRxns > 1) ! number of control variables using OMEGA, not w.   
            ! allocate spate
            allocate ( concIni(nconc, nXpts), concBC(nconc, 2), bioSini(nbioS, nXpts), bioSBC(nbioS, 2) )
            allocate ( cntlIni(ncntlOMG, ncntl(2)) ) 
            allocate ( lowerBnd(ncntl(3)), upperBnd(ncntl(3)) )
            allocate ( tOCnodes(0:ncntl(1)) )
            allocate ( xOCnodes(ncntl(2)) )        
            allocate ( Xpts(nXpts) )
            allocate ( Cc(nconc+nbioS), cLat(nconc+nbioS) )
            allocate ( vSink(nconc+nbioS) )
            allocate ( cntl(0:ncntl(1),ncntl(2),ncntl(3)) )
            allocate ( kIntg(nXpts) )  
            ! read in model parameters
            read(inpUnit,nml=params)
            close(inpUnit)   
         end if
   
         ioerr = abs(ioerr) ! remove negative sign, if present
         ! sum up ioerr across processes and place in ioE in ALL processes.  ioE should be 0 if no errors occured.
         call MPI_ALLREDUCE( ioerr, ioE, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)   
         ! check for error in opening the PRM file by any of the processes.
         if (ioE /= 0) then
            ! an error occured reading the prm file by one or more processes
            if (myRank == 0) then
               write(6,'(/a/)') 'Error opening or reading '//trim(basefilename)//'.inp file. Aborting.'
            end if
            call cleanUp
         end if         
         
         ! Setup uniform points in space for output of solution from BACOLI95
         Xpts = [(dble(i)/dble(nXpts-1)*depth, i=0,nXpts-1)]               
         ! Initialize state and control vectors from data file. Make this read only so all processes can read it.
         open(newunit=icUnit, file=trim(basefilename)//'_IC.dat',action='read',status='old',iostat=ioerr) 
         ! the data file must have the same format as the calling structure below, so this will be problem dependent
         call initArray (icUnit, nconc,    nXpts,    Xpts,     concIni) ! initialize state vector of conentrations
         call initArray (icUnit, nbioS,    nXpts,    Xpts,     bioSini) ! initialize state vector of biological structure
         ! Initial the control vector; however, what is saved in the *_IC.dat file is not what is stored in cntl
         ! because OMEGA is save in the cntl file, NOT w, (lower case omega).  Consequently, nSubRxns set in TypeVar
         ! is used to read in the correct number of cntl variables.
         call initArray (icUnit, ncntlOMG, ncntl(2), xOCnodes, cntlIni) ! initialize control vector
         close(unit=icUnit)
         ! Now convert cntlIni that has OMEGAs to w
         allocate( temp(ncntl(2)) )
         iOMG = 1; iWW = 1
         do i=1,nbioS
            ! Copy the eps value
            cntlIni(iWW,1:ncntl(2)) = cntlIni(iOMG,1:ncntl(2))            
            if (nSubRxns(i) == 1) then
               iOMG = iOMG + 1; iWW = iWW + 1
            else
               temp = 0._dp
               do j=1,nSubRxns(i)-1
                  iOMG = iOMG + 1; iWW = iWW + 1
                  cntlIni(iWW,1:ncntl(2)) = temp + cntlIni(iOMG,1:ncntl(2))
                  temp = temp + cntlIni(iOMG,1:ncntl(2))
               end do
               iOMG = iOMG + 2; iWW = iWW + 1
            end if
         end do  
         deallocate( temp )
         
         ! open up file using MPI for error output from each process.  This is where mpiFHerr gets set in module mpeErrOut.
         ! use MPI_MODE_WRONLY for write only, creat it if necessary, and for sequential access
         ! However, there is no easy way to clear out a file if it already exists, so open with close_on_delete, to get rid of trash, then reopen
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL+MPI_MODE_DELETE_ON_CLOSE, &
                 MPI_INFO_NULL, mpiFHerr, mpiErr)
         call MPI_FILE_CLOSE(mpiFHerr, mpiErr) ! This will delete the file.    
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL, &
                 MPI_INFO_NULL, mpiFHerr, mpiErr)
         ! write a header to the file
         if (myRank == 0) then
            string = 'Errors associated with integration, etc:' 
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
         end if   
         ! Open another file for MPI write for BOBYQA output
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_BOB.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL+MPI_MODE_DELETE_ON_CLOSE, &
                 MPI_INFO_NULL, mpiBOBout, mpiErr)
         call MPI_FILE_CLOSE(mpiBOBout, mpiErr) ! This will delete the file.    
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_BOB.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL, &
                 MPI_INFO_NULL, mpiBOBout, mpiErr)
         ! write a header to the file
         if (myRank == 0) then
            string = 'Variables = "Time" "CPU (min)" "Depth (m)" "Depth Iter" "ojbFcnValue" "Fcn Calls" "errCnt"' 
            call MPI_FILE_WRITE_SHARED(mpiBOBout, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
            string = 'Zone T="BOBYQA output"'
            call MPI_FILE_WRITE_SHARED(mpiBOBout, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
         end if   
         
         if (myRank /=0 ) return ! only the root process will be allowed to write to files (except for error output).
         
         ! open output files and write headers for Tecplot
         iTec = nXpts ! For Tecplot output
         jTec = 0 ! this gets incremented with each save.
         kTec = 1         
         ! Grid points
         open(newunit=gridUnit, file=trim(basefilename)//'_grid.dat',status='unknown')
         string = 'Variables = "Time (d)" "depth (m)"'
         write(gridUnit,'(a)') trim(string)         
         write(gridUnit,'(a)') 'Zone T="Grid Points"'

         open(newunit=OCgridUnit, file=trim(basefilename)//'_OCgrid.dat',status='unknown')
         string = 'Variables = "Time (d)" "depth (m)"'
         write(OCgridUnit,'(a)') trim(string)         
         write(OCgridUnit,'(a)') 'Zone T="OC Grid Points"'
         
         ! Concentration variables (DIRECT ACCESS FILE)
         ! **** Note, for this an the other direct access files, recl must be long enought to contain all the data, eithet character or actual numbers!!!!
         open(newunit=concUnit, file=trim(basefilename)//'_conc.dat',status='unknown', access='direct', recordtype='fixed',recl=240 , form='formatted', &
            carriagecontrol='list') ! Note, this only adds a LF character to the end, but this it standard for linux, and most editors can deal with it.
         ! I have not been able to find how to add CRLF in the Fortran IO (but I did not try the compiler option, fpscomp:general option).
         ! You can add a CRLF to the output, but it is best to add it to the start of a record, such as "write(...) CRLF, number, number" , etc.
         write(concUnit,'(a)',rec=1) 'Variables = "Time (d)" "depth (m)" "Salinity (PSU)" "O2 (mmol/m3)" "DIC (mmol/m3)" "PO4 (mmol/m3)" "SO4 (mmol/m3)" "H2S (mmol/m3)" "C_Phy (mmol/m3)"&
            & "C_GSB (mmol/m3)" "C_L (mmol/m3)" "C_D (mmol/m3)" "P_D (mmol/m3)"'
         write(concUnit,'(a)',rec=2) 'Zone T="Conc vars"'
         write(concUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! Biological structure variables (DIRECT ACCESS FILE)
         open(newunit=bioSUnit, file=trim(basefilename)//'_bioS.dat',status='unknown', access='direct', recordtype='fixed',recl=180 , form='formatted', &
            carriagecontrol='list')
         ! just specifiy the useful names, no need to created numbered ones.
         string = 'Variables = "Time (d)" "depth (m)" "phyS (mmol/m3)" "GSBS (mmol/m3)" "GzS (mmol/m3)" "AGzS (mmol/m3)" "BacS (mmol/m3)" "SRBS (mmol/m3)" "PhS (mmol/m3)" "SOxS (mmol/m3)"' 
         write(bioSUnit,'(a)',rec=1) trim(string)         
         write(bioSUnit,'(a)',rec=2) 'Zone T="BioS vars"'
         write(bioSUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
        
         ! Control variables (DIRECT ACCESS FILE)
         open(newunit=cntlUnit, file=trim(basefilename)//'_cntl.dat',status='unknown', access='direct', recordtype='fixed',recl=600 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "eps-Phy" "omg1-Phy" "omg2-Phy" "eps-GSB" "omg1-GSB" "omg2-GSB" "eps-Gz" "omg1-Gz" "omg2-Gz" "omg3-Gz" "omg4-Gz" "omg5-Gz"&
            & "omg6-Gz" "omg7-Gz" "omg8-Gz" "eps-AGz" "omg1-AGz" "omg2-AGz" "omg3-AGz" "omg4-AGz" "omg5-AGz" "omg6-AGz" "omg7-AGz" "omg8-AGz" "eps-Bac" "omg1-Bac" "omg2-Bac" "omg3-Bac"&
            & "eps-SRB" "omg1-SRB" "omg2-SRB" "omg3-SRB" "eps-Ph" "eps-SOx"' ! just use useful names
         write(cntlUnit,'(a)',rec=1) trim(string)         
         write(cntlUnit,'(a)',rec=2) 'Zone T="cntl vars"'
         write(cntlUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec

         ! reaction variables (DIRECT ACCESS FILE)
         open(newunit=rxnUnit, file=trim(basefilename)//'_rxn.dat',status='unknown', access='direct', recordtype='fixed',recl=620 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "Phy-r1" "Phy-r2" "GSB-r1" "GSB-r2" "Gz-r1" "Gz-r2"&
                & "Gz-r3" "Gz-r4" "Gz-r5" "Gz-r6" "Gz-r7" "Gz-r8" "AGz-r1" "AGz-r2" "AGz-r3"&
                & "AGz-r4" "AGz-r5" "AGz-r6" "AGz-r7" "AGz-r8" "Bac-r1" "Bac-r2" "Bac-r3"&
                & "SRB-r1" "SRB-r2" "SRB-r3" "Ph-r1" "SOx-r1"'
         write(rxnUnit,'(a)',rec=1) trim(string)         
         write(rxnUnit,'(a)',rec=2) 'Zone T="Rxn Rates (mmol/m3/d)"'
         write(rxnUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec

         ! F_K kinetic driver variables (DIRECT ACCESS FILE)
         open(newunit=F_Kunit, file=trim(basefilename)//'_F_K.dat',status='unknown', access='direct', recordtype='fixed',recl=620 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "Phy-r1" "Phy-r2" "GSB-r1" "GSB-r2" "Gz-r1" "Gz-r2"&
                & "Gz-r3" "Gz-r4" "Gz-r5" "Gz-r6" "Gz-r7" "Gz-r8" "AGz-r1" "AGz-r2" "AGz-r3"&
                & "AGz-r4" "AGz-r5" "AGz-r6" "AGz-r7" "AGz-r8" "Bac-r1" "Bac-r2" "Bac-r3"&
                & "SRB-r1" "SRB-r2" "SRB-r3" "Ph-r1" "SOx-r1"'
         write(F_Kunit,'(a)',rec=1) trim(string)         
         write(F_Kunit,'(a)',rec=2) 'Zone T="F_K"'
         write(F_Kunit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! F_T thermodynamic driver variables (DIRECT ACCESS FILE)
         open(newunit=F_Tunit, file=trim(basefilename)//'_F_T.dat',status='unknown', access='direct', recordtype='fixed',recl=620 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "Phy-r1" "Phy-r2" "GSB-r1" "GSB-r2" "Gz-r1" "Gz-r2"&
                & "Gz-r3" "Gz-r4" "Gz-r5" "Gz-r6" "Gz-r7" "Gz-r8" "AGz-r1" "AGz-r2" "AGz-r3"&
                & "AGz-r4" "AGz-r5" "AGz-r6" "AGz-r7" "AGz-r8" "Bac-r1" "Bac-r2" "Bac-r3"&
                & "SRB-r1" "SRB-r2" "SRB-r3" "Ph-r1" "SOx-r1"'
         write(F_Tunit,'(a)',rec=1) trim(string)         
         write(F_Tunit,'(a)',rec=2) 'Zone T="F_T"'
         write(F_Tunit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! reaction free energy variables (DIRECT ACCESS FILE)
         open(newunit=dGrUnit, file=trim(basefilename)//'_dGr.dat',status='unknown', access='direct', recordtype='fixed',recl=620 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "Phy-r1" "Phy-r2" "GSB-r1" "GSB-r2" "Gz-r1" "Gz-r2"&
                & "Gz-r3" "Gz-r4" "Gz-r5" "Gz-r6" "Gz-r7" "Gz-r8" "AGz-r1" "AGz-r2" "AGz-r3"&
                & "AGz-r4" "AGz-r5" "AGz-r6" "AGz-r7" "AGz-r8" "Bac-r1" "Bac-r2" "Bac-r3"&
                & "SRB-r1" "SRB-r2" "SRB-r3" "Ph-r1" "SOx-r1"'
         write(dGrUnit,'(a)',rec=1) trim(string)         
         write(dGrUnit,'(a)',rec=2) 'Zone T="delG (kJ/mol)"'
         write(dGrUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! reaction entropy production (DIRECT ACCESS FILE)
         open(newunit=sigmaRijUnit, file=trim(basefilename)//'_sigmaRij.dat',status='unknown', access='direct', recordtype='fixed',recl=620 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "Phy-r1" "Phy-r2" "GSB-r1" "GSB-r2" "Gz-r1" "Gz-r2"&
                & "Gz-r3" "Gz-r4" "Gz-r5" "Gz-r6" "Gz-r7" "Gz-r8" "AGz-r1" "AGz-r2" "AGz-r3"&
                & "AGz-r4" "AGz-r5" "AGz-r6" "AGz-r7" "AGz-r8" "Bac-r1" "Bac-r2" "Bac-r3"&
                & "SRB-r1" "SRB-r2" "SRB-r3" "Ph-r1" "SOx-r1"'
         write(sigmaRijUnit,'(a)',rec=1) trim(string)         
         write(sigmaRijUnit,'(a)',rec=2) 'Zone T="sigmaRij (J/m/K)"'
         write(sigmaRijUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! sigmaDot and sigma variables (DIRECT ACCESS FILE)
         open(newunit=sigmaUnit, file=trim(basefilename)//'_sigma.dat',status='unknown', access='direct', recordtype='fixed',recl=160 , form='formatted', &
            carriagecontrol='list')
         string = 'Variables = "Time (d)" "depth (m)" "sumSigmaDot (J/m/d/K)" "sumSigma (J/m/K)" "sumSigmaDotRxn (%)" "sumSigmaDotH2O (%)" "sumSigmaDotPart (%)"'
         write(sigmaUnit,'(a)',rec=1) trim(string)         
         write(sigmaUnit,'(a)',rec=2) 'Zone T="Entropy Production"'
         write(sigmaUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
        
         ! store info from interval integration
         open(newunit=summaryUnit, file=trim(basefilename)//'_summary.dat',status='unknown')
         string = 'Variables = "Interval" "Time (d)" "CPU Time (min)" "local iter" "local Residual" "sumSigma Infinite (J/K)" "sumSigma Interval (J/K)" "errCnt"'
         write(summaryUnit,'(a)') trim(string)         
         write(summaryUnit,'(a)') 'Zone T="Summary Data"'
         
         ! Concentration variables (DIRECT ACCESS FILE)
         open(newunit=envUnit, file=trim(basefilename)//'_env.dat',status='unknown', access='direct', recordtype='fixed',recl=110 , form='formatted', &
            carriagecontrol='list')
         write(envUnit,'(a)',rec=1) 'Variables = "Time (d)" "depth (m)" "Light (mmol photons/m2/d)" "PAR (umol photons/m2/s)" "kIntg(z)"'
         write(envUnit,'(a)',rec=2) 'Zone T="Env vars"'
         write(envUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec         
         
         ! ADRsol save.  This is an unformattted file.  See bacoli95_io module above
         ! This is only to output ADRsol.  See read_sol for the input sol file
         open(newunit=solUnit, file=trim(basefilename)//'_tf.sol', form='unformatted')
         
         ! Flux of constituents (DIRECT ACCESS FILE)
         open(newunit=flxUnit, file=trim(basefilename)//'_flx.dat',status='unknown', access='direct', recordtype='fixed',recl=380 , form='formatted', &
            carriagecontrol='list')
         write(flxUnit,'(a)',rec=1) 'Variables = "Time (d)" "depth (m)" "f-Salinity (PSU/d)" "f-O2 (mmol/d)" "f-DIC (mmol/d)" "f-PO4 (mmol/d)"&
            & "f-SO4 (mmol/d)" "f-H2S (mmol/d)" "f-C_Phy (mmol/d)" "f-C_GSB (mmol/d)" "f-C_L (mmol/d)" "f-C_D (mmol/d)" "f-P_D (mmol/d)"&
            & "f-PhyS (mmol/d)" "f-GSBS (mmol/d)" "f-GzS (mmol/d)" "f-AGzS (mmol/d)" "f-BacS (mmol/d)" "f-SRBS (mmol/d)" "f-PhS (mmol/d)" "f-SOxS (mmol/d)"'
         write(flxUnit,'(a)',rec=2) 'Zone T="Total Fluxes"'
         write(flxUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
         
         ! Output for BACOLI95
         open(newunit=BACunit, file=trim(basefilename)//'_BACOLI.dat')
         write(BACunit,'(a)') 'Variables = "Time" "ninit" "remeshings" "init_remeshings" "cold_restarts" "accepted_steps" "BDF_order"'
         write(BACunit,'(a)') 'Zone T="BACOLI95 Output:'
         
         return
      end subroutine initialize
      
      subroutine cleanUp
         ! close units where necessary and deallocate
         use mpiErrOut
         integer mpierr
         deallocate ( concIni, concBC, bioSini, bioSBC )
         deallocate ( lowerBnd, upperBnd )         
         deallocate ( tOCnodes ) 
         deallocate ( xOCnodes )        
         deallocate ( Xpts, Cc, cLat )     
         deallocate ( vSink )       
         deallocate ( cntl )
         deallocate ( kIntg )
         if (myRank == 0) then 
            ! Units are only opened in the root process
            close(inpUnit    )
            close(concUnit   )
            close(bioSUnit   )
            close(cntlUnit   )
            close(rxnUnit    )
            close(F_Kunit    )
            close(F_Tunit    )
            close(dGrUnit    )
            close(sigmaRijUnit)   
            close(sigmaUnit  )
            close(summaryUnit)
            close(gridUnit   )
            close(OCgridUnit )
            close(envUnit    )
            close(solUnit    )
            close(flxUnit    )
            close(BACunit    )
         end if
         call MPI_FILE_CLOSE(mpiFHerr, mpiErr)
         call MPI_FILE_CLOSE(mpiBOBout, mpiErr) ! added on 1-Mar-2018             
         call MPI_FINALIZE(mpierr)
         stop         
      end subroutine cleanUp
      
      subroutine  kIntegrate (nC, Cout)
         ! The state solution to the PDEs at the previous time point is used to estimate the integral of the 
         ! light absorption term to approximate the light profile and its derivative since BACOLI is unable to handle PDEs mixed with ODEs
         implicit none
         integer , intent(in):: nC
         real(dp), intent(in):: Cout(nC,nXpts,2) ! Concentration of all state variables from BACOLI at Xpts locations. Note, derivative part is not used here
         ! local declarations
         integer i
         
         ! Integrate light asborbing terms and save cumulatively in kIntg
         kIntg(1) = 0._dp ! integrated value of the light absorption at Xpts.  It includes the negative sign.
         ! only consider detritus as carbon to intercept light (don't double count P_D)
         Associate ( C_Phy => Cout(7,1:nXpts,1), C_GSB => Cout(8,1:nXpts,1), C_D => Cout(10,1:nXpts,1), PhyS => Cout(12,1:nXpts,1), GSBS => Cout(13,1:nXpts,1), &
                     GzS => Cout(14,1:nXpts,1), AGzS => Cout(15,1:nXpts,1), BacS => Cout(16,1:nXpts,1), SRBS => Cout(17,1:nXpts,1), PhS => Cout(18,1:nXpts,1), SOxS => Cout(19,1:nXpts,1) )
            do i=2, nXpts
               kIntg(i) = kIntg(i-1) - ( ( k_w + k_p*(C_Phy(i-1)+C_GSB(i-1)+C_D(i-1)+PhyS(i-1)+GSBS(i-1)+GzS(i-1)+AGzS(i-1)+BacS(i-1)+SRBS(i-1)+PhS(i-1)+SOxS(i-1)) &
                        + k_w + k_p*(C_Phy(i)+C_GSB(i)+C_D(i)+PhyS(i)+GSBS(i)+GzS(i)+AGzS(i)+BacS(i)+SRBS(i)+PhS(i)+SOxS(i)) )/2._dp * (Xpts(i)-Xpts(i-1)) )
            end do
         end associate
         return
      end subroutine kIntegrate   
      
      real(dp) function kIntg_x (x)
         ! This interpolate kIntg at location x
         implicit none
         real(dp) x
         ! local declarations         
         integer nlow
         real(dp) delt
         
         call interp1D (x, nXpts, Xpts, nlow, delt)
         kIntg_x = kIntg(nlow) + delt*(kIntg(nlow+1) - kIntg(nlow))
         return
      end function kIntg_x            
   end module parameters
      
   module functions
      ! Note, functions are placed in this module so that real(dp) can be used instead of real(8), etc.
      ! It's not really necessary...
      use realPrec
      implicit none      
   contains
      real(dp) function I0(t)
         ! This fuction returns the light at the surface (mmol photons /m^2 /d (Not micromoles)) at time t (d)
         use parameters, only: I0max, dLat
         real(dp) t
         real(dp) tJulian! since t is just linear time, it needs to be converted to a julian day, but don't bother with leap years.
         tJulian = mod(t,365._dp)
         call parSur (tJulian, dLat, I0max, I0)
         
         return
      end function I0
      
      real(dp) function weightS (t)
         ! weighting (or discouting) function on sigma
         use parameters, only: wInfinity, deltInfinity, useSigmaW, tOCnodes
         implicit none
         real(dp) t ! current time
         ! Local declarations
         real(dp) kW
         kW = -log(wInfinity)/deltInfinity
         weightS = exp(-kW*(t-tOCnodes(0)))
         if (.not. useSigmaW) weightS = 1.0_dp ! When solving over interval, don't weight
         return
      end function weightS      
      
      real(dp) function F_Thermo (delG, ne, Tk)
         ! This function calculates the thermodynamic driver using LaRowe2012
         ! Function returns a unitless value between 0 and 1
         use parameters, only: delPsi  ! The electric potential accross the cell wall (V).  Note, this is sigmodial 
                                       ! function, so F_Thermo is not zero if delG is less than it.  See LaRowe2012
         real(dp) delG  ! free energy of reaction (kJ/mol reaction)
         real(dp) ne    ! mole of electrons exchanged per mole of reaction extent
         real(dp) Tk    ! Temperature (K)
         ! Local declarations
         real(dp), parameter:: RkJ  = 8.3144598d-3   ! gas constant (kJ/(g-mol K) or J/(g-mmol K))
         real(dp), parameter:: F = 96485.3329 ! Faraday constant (C/mol-e)
         
         F_Thermo = 1.0_dp/( 1.0_dp + Exp( (delG/ne + F*delPsi/1000._dp)/(RkJ*Tk) ) )  ! Divide F*delPsi to get kJ/mol instead of J/mol
      end function  F_Thermo   
      
      real(dp) function stepUp(x, xm, sig)
         ! This routine generates a continuous sigmodial step up function from 0 to 1 around xm, with a gradient described by sig
         implicit none
         real(dp) x   ! value of x where function is evaluated
         real(dp) xm  ! value of x where setUp is 0.5
         real(dp) sig ! how steap the exponential function is around xm
         
         stepUp = 1._dp/( 1._dp + exp(-sig*(x - xm)) )
         return
      end function stepUp      
      
      real(dp) function fZero(absZero, x)
         ! This routine is used to prevent state variables from going negative or to zero, but in a continuous manner
         ! That is, while x >= 2 xm, the function returns x, but as x decreases below 2 xm, fZero also returns an number: xm < fZero <= 2 xm
         implicit none
         real(dp) absZero  ! Value that x must stay above or at
         real(dp) x        ! value of x that could be zero or less than zero that needs to be prevented from that

         if (x > 2._dp*absZero) then
            fZero = x
            return
         end if
         if (x < 0._dp) then
            fZero = absZero
            return
         end if
         ! x is between 0 and 2 absZero, so use quadratic for transition
         fZero = x**2/(4._dp*absZero) + absZero
         return
      end function fZero      
   end module functions

   program MEPOCphototroph_1D
      ! This program is a Optimal Control 1D model to investigate Siders Pond
      ! This is the MPI version using BOBYQA
      use realPrec
      use parameters
      use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown    
      use mpi
      use mpiErrOut, only: mpiFHerr ! Unit for writing errors to for bacoli (see bacoli_errorOut.f)
      use omp_lib
      implicit none
      real(dp) tInt0, tIntf
      integer i, ifail, ifailTot, nIntervals, iostat, noProc, mpierr, maxThreads, noThreads
      real(dp), allocatable:: sumSigmaInfinite(:), sumSigmaInterval(:) 
      real(dp) sumSigmaTotInf, sumSigmaTotFin
      integer IUNITA(5) ! needed by SLATEC error handling routine.  See call to XSETUA below
      integer nameLen
      character (len=MPI_MAX_PROCESSOR_NAME) nodeName
      character(len=8) timeStr
      character(len=9) dateStr

      ! Initialize MPI
      call MPI_INIT( mpierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpierr) ! get rank of this process in world

      if (myRank == 0) write(*,'(/a/)') 'Program with BOBYQA Version 1.5.9, 31-Dec-2017'
      call initialize
      ! routines in bacoli-aux use SLATEC error handling. Eventhough a common block is used, it is only written to once
      ! here by XSETUA, so does not pose a real problem for threads.
      IUNITA(1) = mpiFHerr
      call XSETUA (IUNITA, 1)
      
      ! Get the number of processes running and the maxium number of threads
      call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)
      maxThreads = OMP_GET_MAX_THREADS()
      call omp_set_num_threads(min(maxThreads, ompThreads)) ! Don't use more threads than available
      !$OMP parallel
      noThreads = OMP_GET_NUM_THREADS() ! Outside of omp parallel regions, this returns 1
      !$OMP end parallel
      if (myRank == 0) write(*,'(a,i4,a)') 'Running program with ',noProc, ' processes'     
      call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      call MPI_GET_PROCESSOR_NAME(nodeName, nameLen, mpiErr)
      do i=0,noProc-1
         if (i == myRank) then
            write(*,'(a,i4,2a,2(a,i4))') ' Process ', myRank, ' on node "', trim(nodeName), '" has ', noThreads, ' OpenMP threads out of: ', maxThreads
            write(*,'(a,i12,/)')  ' The thread stack size is set to: ', KMP_GET_STACKSIZE_S()   
         end if
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      end do      
      
      WTime_0 = MPI_Wtime() ! Start the time clock
      ! setup initial control variables. Matrix takes the form cntl(nt, nx, nv), but cntl(0,1:nx,1:nv) stores solution at 0 from previous interval.
      do i=0,ncntl(1)
         cntl(i,1:ncntl(2),1:ncntl(3)) = transpose( cntlIni(1:ncntl(3),1:ncntl(2)) ) 
      end do      
      deallocate ( cntlIni ) ! this is no longer needed now that they have been copied to cntl
      allocate ( sumSigmaInfinite(ncntl(2)), sumSigmaInterval(ncntl(2)) ) 

      tOCnodes(0) = tStart ! set the initial time
      call spaceNodes () ! This set the time points in tOCnodes to approximate the control variables over deltInfinity
      solGood = .false. ! whether ADRsol has a valid solution.   
      call getADRsol (iostat)  ! read in ADRsol if readADRsol is true. 
      if (iostat /= 0) then
         ! Since I may not want to mess up file, just stop on io error
         if (myRank == 0) write(*,'(a,i)') 'Error reading ADRsol from *.sol file. IOSTAT = ', iostat
         call cleanUp
         stop
      end if     
      if (solGood .and. myRank == 0) write(*,'(a, g13.5,/)') 'Using stored ADRsol solution with t0 = ', ADRsol%t0
      ! Save iniatial solution either based on Cinit or ADRsol
      call saveInitialSolution (tStart)
      ! search for solutions over all intervals, but if PDEonly true, then don't change cntl and just run
      nIntervals = ceiling((tFinish - tStart)/deltInterval) ! Number of intervals needed to reach tEnd (solution may go beyond tEnd)
      Do interval=1, nIntervals
         call spaceNodes () ! This set the time points in tOCnodes to approximate the control variables over deltInfinity
         sumSigmaTotInf = 0._dp ! just zero out in case PDEonly set to true (doesn't really matter though)
         errCnt = 0 ! set the number of errors to zero for interval, for this process
         ! Main optimization loop occurs within this if construct
         if (.not. PDEonly) then
            saveSoln_BACOLI = .false.   ! don't save ODE/PDE solution while searching for maximum.
            useSigmaW = .true. ! discount sigma solution with weighting
            call findMaxEntropy (sumSigmaInfinite, sumSigmaTotInf, ifail) 
         end if
         ! Use cntl solution and integrate over the update interval, saving solution
         tInt0 = tOCnodes(0)
         tIntf = tInt0 + deltInterval ! Integrate only up to end of interval
         saveSoln_BACOLI = .true.    ! save solution over interval
         useSigmaW = .false. ! don't discount sigma solution
         call saveOCgrid (tIntf) ! save the node points for the current finite interval
         call integrateState(tInt0, tIntf, sumSigmaInterval, sumSigmaTotFin, ifail)
         ! If any process fails, then stop all proceses.
         call MPI_ALLREDUCE( ifail, ifailTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum failures for all processes.           
         if (ifailTot /= 0) then    
            if ( myRank == 0 ) then
               write(*,'(a)') 'Error::integration over interval failed for some process during solution storage!'
               write(*,'(2(a,f6.1))') ' Got to time: ', tIntf, ' instead of tIntf: ', tInt0 + deltInterval
               write(*,'(a)') ' Stopping'
            end if            
            call cleanUp
         end if
         if (myRank == 0 .and. solGood .and. ifail /= 1) then
            rewind (solUnit) ! only save the last good solution.
            call write_sol(solUnit, ADRsol, iostat) ! saves solution so it can be recovered
         end if   
         call MPI_ALLREDUCE( errCnt, errCntTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum errors for all processes.  
         call summarySave(interval, tIntf, sumSigmaTotInf, sumSigmaTotFin)
         tOCnodes(0) = tIntf ! copy time of end of interval to the start of the next. In theory, this should work for all
         ! copy control variables at end of interval to IC for next interval
         if (.not. PDEonly) cntl(0,:,:) = cntl(ncntl(1),:,:)  ! all other values are not changed, just use last solution for guess.   
      end do
      WTime_f = MPI_Wtime() ! end timer      
      call time(timeStr); call date(dateStr)
      if (myRank == 0) write(*,'(a,f9.2,a)') 'Finished::Execution clock time (min): ', (WTime_f-WTime_0)/60._dp, ' at '//timeStr//' on '//dateStr
      
      ! clean things up locally
      call bacoli95_sol_teardown(ADRsol)      
      deallocate ( sumSigmaInfinite, sumSigmaInterval )
      ! Deallocate global vars, close units, and stop
      call cleanUp
   end program MEPOCphototroph_1D
   
   subroutine initArray (iUnit, nVar, nX, X, vec)
      ! this routine is used to initialize an array based on data read in
      ! This has been update to treat any record begining with ! as a comment.  
      ! However, comments CANNOT be at ends of input lines, NOR can they be inserted
      ! within the data block as given below.  That is, comments are only allows before a 
      ! a data block
      use realPrec
      implicit none
      integer, intent(in)  :: iUnit         ! The unit that has already been opened for read
      integer, intent(in)  :: nVar          ! number of variables in vec, which does not include X (or t that is dropped)      
      integer, intent(in)  :: nX            ! size of X and vec
      real(dp), intent(in) :: X(nX)         ! spatial variable to interpolate over
      real(dp), intent(inout):: vec(nVar, nX) ! vector or array to store dependent variable in
      ! The general data format of a file alread opend in unit iUnit is
      !   nDx, nDv          This is the number of variables in the data file and nDx <= nX and nDv <= nY
      !   t, x, var1, var2, ... nVar     Currently, the t value is just dropped
      !   ...
      !   nX records
      ! Local declarations
      integer i, j, nDx, nDv, nlow
      integer cmt
      real(dp) t, delx, junk
      real(dp), allocatable:: dX(:), dataVX(:,:)
      character(len=300) rec ! This means the maximum record length is 300, which seem big enough.
      
      do
         ! Read through comment lines until nDx, nDy is read
         read(iUnit,'(a)') rec
         if (len_trim(rec) /= 0 .and. scan(rec,'!') == 0) exit
      end do
      read(rec,*) nDx, nDv
      if (nDx < 2 .or. nDv /= nVar) then
         write(*,*) 'Error reading data file for initialization: nDx < 2 or nDv /= nVar'
         stop
      end if     
      allocate ( dX(nDx), dataVX(nDv, nDx) )  ! this is allocate in transpose to how data is in the file
      do i=1,nDx
         read(iUnit,*) t, dX(i), (dataVX(j,i), j=1,nDv)
      end do
      ! using this data, interpolate to vec using x
      do i=1,nX
         call interp1D (X(i), nDx, dX, nlow, delx)
         vec(:,i) = dataVX(:,nlow) + delx*( dataVX(:,nlow+1) - dataVX(:,nlow) )
      end do
      deallocate ( dX, dataVX )
      return
   end subroutine initArray
   
   subroutine findMaxEntropy (sumSigmaX, sumSigmaTot, ifail)
      ! this routine finds the maximum entropy production over an interval, but using
      ! local optimization for the spatial dimension, since phototrophs populations will
      ! only match reality if entropy is max locally.
      ! This version uses BOBYQA and parallelizes over depths with MPI
      use realPrec
      use parameters
      use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown, bacoli95_init  
      use mpi
      use interfaces
      implicit none
      
      ! Dummy variables
      real(dp), intent(out)  :: sumSigmaX(ncntl(2)) ! total internal entropy production over interval at xOCnodes.
      real(dp), intent(out)  :: sumSigmaTot         ! total entropy produced over domain.
      integer, intent(out)   :: ifail               ! Set to 1 if error occured (right now, just in integration)
      
      ! other declarations
      integer nU ! number of unknowns
      integer i, j, iX, iter, mpierr, noProc
      real(dp) sumSigmaXold(ncntl(2)) ! used to hold the previous iteration of sumSigmaX
      character*8 timeStr
      character(len=9) dateStr
      character(len=80) fmtstr, string
      real(dp) cntlPrev(ncntl(1), ncntl(2), ncntl(3)), endTi, startTi
      real(dp), allocatable  :: lowerU(:), upperU(:), w_bob(:)
      
      ifail = 0 ! set error flag (=1 if error occurs)
      nU = ncntl(1)*ncntl(3)  ! at each depth in xOCndoes, there are this many variables to be adjusted
      ! Allocate dynamic arrays 
      allocate ( lowerU(nU), upperU(nU) )
      allocate ( w_bob((npt_bobyqa+5)*(npt_bobyqa+nU)+3*nU*(nU+5)/2) )
      
      ! Conduct local entropy maximization at each xOCnode
      ! Allocate dynamic arrays.
      lowerU(1:nU)  = [ ((lowerBnd(j),i=1,ncntl(1)),j=1,ncntl(3)) ]  ! all tOCnodes and xOCnodes have the same boundaries on cntl vars.
      upperU(1:nU)  = [ ((upperBnd(j),i=1,ncntl(1)),j=1,ncntl(3)) ] 
      
      ! Get the number of processes running
      call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)

      ! Begin loop to locally maximize entropy production at each xOCnode location.
      cntlPrev = cntl(1:ncntl(1), 1:ncntl(2), 1:ncntl(3)) ! used to measure changes in cntl between iterations
      ! Begin muliple refinements of solution
      do iter = 1, maxLocalIter
         call time(timeStr); call date(dateStr)
         if (myRank == 0) write(*,'(a,i3,a)') 'Starting depth interation ', iter, ' of local optimization at time: '//dateStr//' '//timeStr
         ! Begin MPI parallel loop.  Although master and slave might be better, for now just do a simple loop in which
         ! each process only executes part of the loop.
         ! In this case, as each depth is optimized, the new control variables at the depth are used in the next.
         ! In the current simple do parallization, the solution will only be "reproducible" if the same number of processes are run
         ! For instance, if just one process is used, then cntl(:,i,:) will be updated from top to bottom.  If there are as many
         ! processors as depths, then cntl remains unchanged execpt for the depth being optimized.
         do iX =  myRank+1, ncntl(2), noProc
            iNode = iX ! this gets passed to calfun, the objective function  
            fcnCalls = 0
            startTi = MPI_WTime()
            ! Call BOBYQA
            call BOBYQA (nU ,npt_bobyqa , cntl(1:ncntl(1), iX, 1:ncntl(3)), lowerU, upperU, rhobeg, rhoend, iprint, maxfun, w_bob)            
            endTi = MPI_WTime()            
            ! Get value of entropy at maximum for current node
            call entropy(sumSigmaX, SumSigmaTot, ifail)
            ! Store output from BOBYQA
            fmtstr = '(f10.4,1x,f7.2,1x,f5.1,1x,i3,1x,g13.5,1x,i6,1x,i6)'
            write(string,fmtstr) tOCnodes(ncntl(1)), (endTi-startTi)/60.0, xOCnodes(iNode), iter, sumSigmaX(iX), fcnCalls, errCnt   
            call MPI_FILE_WRITE_SHARED(mpiBOBout, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
            call time(timeStr)
            ! Write to monitor
            write(*,'(a,f5.1,a,i3,a,a,G13.5,a,i5)') 'Search at depth: ', xOCnodes(iX), ' on iteration: ', iter, &
               ' time: '//timeStr, ' sumSigma: ', sumSigmaX(iX), ' errCnt: ', errCnt
         end do
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr) ! this is probably not necessary, but wait for all processes to finish their jobs
         ! broadcast the solutions found by this process to all others so that cntl gets updated for the next iteration or interval
         ! This process has solved for those depths starting a myRank+1 with a stride of noProc
         do i = 1, ncntl(2)     
            ! In MPI-BCAST both the sending and receiving processes must specify the same root.
            call MPI_BCAST(cntl(1:ncntl(1), i, 1:ncntl(3)), nU, MPI_DOUBLE_PRECISION, mod(i-1, noProc), MPI_COMM_WORLD, mpierr)
         end do
         localIter = iter ! Number of times the local solution has been attempted
         ! Determine if iterations have converged to local maximums based on changes in control variables.
         localRes = sum(abs(cntl(1:ncntl(1), 1:ncntl(2), 1:ncntl(3)) - cntlPrev(1:ncntl(1), 1:ncntl(2), 1:ncntl(3))))/product(ncntl)
         if (localRes <= maxLocalRes) exit ! local maximum iterations have converged (but there may be more tests that
         ! could be done, such as comparing changes in cntl between iterations or max(iX) localy vs at end of iteration. 
         cntlPrev = cntl(1:ncntl(1), 1:ncntl(2), 1:ncntl(3)) ! copy to previous 
      end do
      ! Return SumSigmaTot for the final value of cntl, where cntl should be the same for all processes. 
      call entropy(sumSigmaX, SumSigmaTot, ifail)
      deallocate ( lowerU, upperU )
      deallocate ( w_bob )
      return
   end subroutine findMaxEntropy
   
   subroutine calfun (nU, x, f)
      ! this routine is called by BOBYQA to find the function value, in this case entropy at a given level.
      use realPrec
      use parameters
      implicit none
      ! Dummy variables
      integer, intent(in)  :: nU    ! number of control variables to be optimized																		
      real(dp), intent(in) :: x(nU) ! vector of control values at a given depth
      real(dp), intent(out):: f     ! Entropy production at current depth
      ! Local declarations
      real(dp) sumSigmaX(ncntl(2))  ! entropy production is calculated at each spatial node, but only need the current one
      real(dp) SumSigmaTot          ! entropy production integrated over depth, but not used here
      integer ifail                 ! if there was a problem with integration.
      integer i, j
      
      ! Update cntl with values passed by BOBYQA, which updates cntl globally via parameters module
      do j=1,ncntl(3)
         do i=1,ncntl(1)
            cntl(i, iNode, j) = x(i+ncntl(1)*(j-1))
         end do
      end do
      fcnCalls = fcnCalls + 1
      call entropy(sumSigmaX, SumSigmaTot, ifail)      
      f = -sumSigmaX(iNode)  ! BOBYQA finds minimums, so take negative.
      if (ifail /= 0) f = huge(f) ! BOBYQA does not handle failures, so just set to huge
      return
   end subroutine calfun    
   
   subroutine entropy(sumSigmaX, sumSigmaTot, ifail) 
      ! This routine calculates the total entropy production over the interval
      ! Note, this routine is embeded in an openMP construct started in entropyWgradient
      use realPrec
      use parameters
      use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown, bacoli95_init            
   
      implicit none
      real(dp), intent(out):: sumSigmaX(ncntl(2))  ! sum of internal entropy production at the xOCnodes
      real(dp), intent(out):: sumSigmaTot          ! Total internal entropy production (integration of sumSigmaX)
      integer,  intent(out):: ifail ! if integration failed ifail > 0
      ! local declarations
      integer i
      real(dp) t0, tf
      type(bacoli95_sol) ADRsolLoc ! local copy
      logical solGoodLoc
 
      t0 = tOCnodes(0); tf = tOCnodes(ncntl(1)) ! This can change, so do not modify tOCnodes
      solGoodLoc = solGood 
      ! Need to setup a copy ADRsol, to restore the global state because call to Entropy is not being used to advance the PDE integration.
      if ( solGood ) then
         call bacoli95_init(ADRsolLoc, nconc+nbioS, Xpts(1:nXpts)/Xc, atol=atol1_bacol, rtol=rtol1_bacol, &
            nint_max=nXpts_max, kcol=kcol_bacol) ! this allocates space 
         call copySol (ADRsol, ADRsolLoc) ! copy values not pointers
      end if      
      call integrateState(t0, tf, sumSigmaX, sumSigmaTot, ifail)
      solGood = solGoodLoc ! return solGood to it's previous state
      if ( solGood ) then
         ! return ADRsol to its previous state
         call copySol (ADRsolLoc, ADRsol) ! copy values not pointers
         call bacoli95_sol_teardown(ADRsolLoc) 
      end if            
      return
   end subroutine entropy
   
   subroutine integrateState(t0, tf, sumSigmaOC, sumSigmaTot, ifail)
      ! This routine integrates the state PDE equations using BACOLI95
      ! Most parameters are passed via module parameters.
      ! Also, the initial conditions, which are given in the routine 
      ! In this version, sigma is not added to the set of PDEs, but it is kept in the solution vector.
      use realPrec
      use parameters
      use functions, only: weightS
      use mpi
      use mpiErrOut
      use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals, &
                              bacoli95_sol, bacoli95_sol_teardown    
      use typeVars
      implicit none
      real(dp), intent(in)             :: t0 ! Start time of integration
      real(dp), intent(inout)          :: tf ! Stop time of integration.  If a error occurs this will be set to where the integration got to.
      real(dp), intent(out)            :: sumSigmaOC(ncntl(2))  ! value of sumSigma at tf at the xOCnodes locations (J/K/m)
      real(dp), intent(out)            :: sumSigmaTot  ! value of sumSigma at tf integrated over xPts locations (J/K)
      integer,  intent(out)            :: ifail       ! if the integration failed (=1); otherwise 0

      ! Local declarations needed for BACOLI95
      integer nC ! number of PDEs
      real(dp) Cout(nconc+nbioS,nXpts,2) ! array to store solution in and it's first spatial derivative.
      character(len=80) string
      external PDEs, bndxa, bndxb, Cinit
      integer nTpts ! number of time points outputs between t0 and tf

      ! Other local declarations   
      integer i, j, nlow, ndcdx, mpiErr
      real(dp) t, delx
      real(dp) sumSigmaX(nXpts), tPrev, sumSigmaDot(3), sumSigmaDotPrev(nXpts), cntlTX(ncntl(3))
      type(rxnNet) rNet ! stores reactions, free energies, etc calculated in rxnProperties, and now also holds control variables at t and x.
   
      ifail = 0 ! set flag to 1 if integration error occurs.
      sumSigmaTot = 0._dp ! initialize here in case integration fails      
      ! Make sure ADRsol makes sense if using
      if (solGood) then
         if (ADRsol%t0 /= t0) then
            write(string,'(a)') 'Warning:: ADRsol time does not match t0 in integrateState'
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
         end if
      end if
      
      ! Space needed by BACOLI95
      ! Set size of PDEs, NOT including sigma integration
      nC = nconc+nbioS     
      ! setup kIntg with either Cinit or from ADRsol
      if (solGood) then
         ! use ADRsol to get Cout
         call bacoli95_vals(ADRsol, Xpts/Xc, Cout(1:nC,1:nXpts,1), nderiv=0) ! get solution at Xpts, but don't need derivatives
         do i=1,nXpts
            Cout(1:nC,i,1) = Cout(1:nC,i,1)*Cc(1:nC) ! dimensionalize output.  Note, derivatives not needed here for kIntegrate 
         end do         
      else
         ! Use Cint to set Cout
         do i=1,nXpts
            call Cinit(Xpts(i)/Xc, Cout(1:nC,i,1), nC)  
            Cout(1:nC,i,1) = Cout(1:nC,i,1)*Cc(1:nC) ! dimensionalize
         end do 
      end if      
      ! Use Cout to inialize kIntg
      ! But first insure state variables are non negative
      where ( Cout(1:nC,1:nXpts,1) < 0._dp ) Cout(1:nC,1:nXpts,1) = 0._dp
      call kIntegrate (nC, Cout) ! initialize kIntg

      ! BACOLI Initialization if solGood is false
      ! The OC spatial grid and the PDE spatial grid should have the same start and end boundary points, but 
      ! intermediate points in x will be different between the OC grid and the PDE grid.  Same holds for time too.
      ! not sure where nint gets set, but I'm guessing it is a default
      ! Only using scaler values for atol and rtol
      if (.not. solGood) then
         call bacoli95_sol_teardown(ADRsol) ! deallocate ADRsol if solGood is false
         call bacoli95_init(ADRsol, nC, Xpts(1:nXpts)/Xc, tstart=t0, atol=atol1_bacol, rtol=rtol1_bacol, &
            nint_max=nXpts_max, kcol=kcol_bacol)
      end if
      
      ! Get summed entropy production at t0 at the nXpts points
      tPrev = t0 ! used for trapizoid rule for getting internal entropy producton.
      do j=1,nXpts
         ! Calculate sumSigmaDot at depth Xpts(j) at t0
         call interp2D (ncntl(1)+1, ncntl(2), ncntl(3), tOCnodes, xOCnodes, cntl, t0, Xpts(j), cntlTX) ! get control vars at (t0,x)
         call rxnProperties(t0, Xpts(j), Cout(1:nconc,j,1), Cout(nconc+1:nconc+nBioS,j,1), ncntl(3), cntlTX, rNet, sumSigmaDot)
         sumSigmaDotPrev(j) = sum(sumSigmaDot)
      end do
      sumSigmaX = 0._dp ! zero out vector used in trapezoid rule.  Starts at zero at start of interval
      
      ! Begin solution integration, storing (if saveSoln_Bim is true) points at equal time intervals
      ! it may be worth exploring if the PDEs are solved faster or differently if intermediate time points are requested  
      ! versus just solving to t
      nTpts = ceiling((tf - t0)/deltaT) ! round up, so the real spacing will be smaller than deltaT
      do i=1, nTpts
         t = t0 + real(i)*(tf-t0)/real(nTpts)
         ! Compute solution at current t
         call bacoli95(ADRsol, t, PDEs, bndxa, bndxb, Cinit)
         if (ADRsol%idid <= 0) then
            ! If this happens very often, then better error handling will need to be placed here.
            errCnt = errCnt + 1
            write(string, '(a,i4,a,f7.2,a,i4)') 'Error in BACOLI::idid = ', ADRsol%idid, ' at time: ', ADRsol%t0, ' for rank: ', myRank
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
            tf = ADRsol%t0 ! set tf to where the solution got to.
            ifail = 1
            return
         end if
         ndcdx = 0 ! no deriviatives
         if (saveSoln_BACOLI) ndcdx = 1 !  need derivatives for solution saving
         call bacoli95_vals(ADRsol, Xpts/Xc, Cout(1:nC,1:nXpts,1:ndcdx+1), nderiv=ndcdx) ! get solution and derivative, if needed, at Xpts
         do j=1,nXpts
            Cout(1:nC,j,1) = Cout(1:nC,j,1)*Cc(1:nC) ! dimensionalize output  
            if (saveSoln_BACOLI) Cout(1:nC,j,2) = Cout(1:nC,j,2)*Cc(1:nC)/xC ! dimensionalize output, derivatives are needed below in solution save
            ! Calculate sumSigmaDot at depth Xpts(j)
            call interp2D (ncntl(1)+1, ncntl(2), ncntl(3), tOCnodes, xOCnodes, cntl, t, Xpts(j), cntlTX) ! get control vars at (t,x)
            call rxnProperties(t, Xpts(j), Cout(1:nconc,j,1), Cout(nconc+1:nconc+nBioS,j,1), ncntl(3), cntlTX, rNet, sumSigmaDot)
            ! Add to integration term, apply future discouting if turned on.
            sumSigmaX(j) = sumSigmaX(j) + weightS(t)*(sum(sumSigmaDot) + sumSigmaDotPrev(j))*(t - tPrev)/2._dp
            sumSigmaDotPrev(j) = sum(sumSigmaDot)
         end do
         tPrev = t         
         
         ! Get solution at t for nXpts values given in Xpts and save it if saveSoln_BACOLI is true
         if (saveSoln_BACOLI) then
            call saveSolution (t, nC, Cout, sumSigmaX)
            call saveGrid (ADRsol)
         end if
         ! setup kIntg for the next time point at the nXpts points
         where ( Cout(1:nC,1:nXpts,1) < 0._dp ) Cout(1:nC,1:nXpts,1) = 0._dp
         call kIntegrate (nC, Cout) ! initialize kIntg                  
      end do
      solGood = .true.  ! was able to integrate through interval, so ADRsol can be used to start next interval.    
      
      ! calculate the total entropy production integrated over depth (J/k)      
      do i=1,nXpts-1
         sumSigmaTot = sumSigmaTot + (sumSigmaX(i+1) + sumSigmaX(i))*(Xpts(i+1) - Xpts(i))/2.0_dp ! dimensionalize
      end do

      ! get value of sumSigmaOC at tf at the xOCnodes locations by interpolation
      do i=1,ncntl(2)   
         call interp1D (xOCnodes(i), nXpts, Xpts, nlow, delx)
         sumSigmaOC(i) = sumSigmaX(nlow) + delx*(sumSigmaX(nlow+1) - sumSigmaX(nlow))
      end do
      return
   end subroutine integrateState
   
   subroutine saveInitialSolution (t0)
      ! This routine saves the initial solution at t0
      use realPrec
      use parameters
      use bacoli95_mod, only: bacoli95_vals
      implicit none
      real(dp), intent(in):: t0 ! Time as start of run for entire solution
      ! local declarations
      integer nC, i
      real(dp) Cout(nconc+nbioS,nXpts,2) ! array to store solution in.
      real(dp) sumSigmaX(nXpts) ! need vector to pass, but will be all zeros at initial time
      
      if (myRank /= 0) return ! only done by the master process
      ! allocate Cout
      nC = nconc+nbioS ! number of PDEs 
      if (solGood) then
         ! use ADRsol to populate Cout
         if (ADRsol%t0 /= t0) then 
            write(*,'(a,g15.5)') 'WARNING t0 in ADRsol does not match t0 in *.inp file. ADRsol%t0 = ', ADRsol%t0 
            pause 'Hit CR to continue or CNTL-C to terminate' ! this happens almost instantly at program start, so ask.
         end if         
         call bacoli95_vals(ADRsol, Xpts/Xc, Cout, nderiv=1) ! get solution and derivatives at Xpts
         do i=1,nXpts
            Cout(1:nC,i,1) = Cout(1:nC,i,1)*Cc(1:nC) ! dimensionalize output 
            Cout(1:nC,i,2) = Cout(1:nC,i,2)*Cc(1:nC)/xC ! dimensionalize output 
         end do       
      else         
         ! populate Cout using Cinit for BALCOLI instead.  Note, this will only occur when the program is first run, AND ADRsol is not read in.
         do i=1,nXpts
            call Cinit(Xpts(i)/Xc, Cout(1:nC,i,1), nC)
            Cout(1:nC,i,1) = Cout(1:nC,i,1)*Cc(1:nC) ! dimensionalize
            Cout(1:nC,i,2) = 0._dp  ! While I could try to calculate derivatives, it's not that important here, set to zero.
         end do  
      end if
      where ( Cout(1:nC,1:nXpts,1) < 0._dp ) Cout(1:nC,1:nXpts,1) = 0._dp      
      call kIntegrate (nC, Cout) ! initialize kIntg for light calculations
      ! Unless I want to write a routine the estimate the derivative of C wrt x, just set those to zero
      ! Save solution
      sumSigmaX = 0._dp
      call saveSolution (t0, nC, Cout, sumSigmaX)      
      return
   end subroutine saveInitialSolution
   
   subroutine saveSolution (t, nC, Cout, sumSigmaX)
      ! save the solution output
      use realPrec
      use parameters
      use functions, only: I0
      use realPrec
      use typeVars
      use functions, only: stepUp
      implicit none
      real(dp), intent(in)::  t                  ! time (d)
      integer,  intent(in)::  nC                 ! row size dimenstions of Cout (number of PDEs)
      real(dp), intent(in)::  Cout(nC, nXpts, 2) ! concentration of constituent variables And derivatives from BACOLI at time t
      real(dp), intent(in)::  sumSigmaX(nXpts)   ! internal entropy production at time t over the nXpts points
      
      ! Local declarations      
      real(dp) sumSigmaDot(3), sumSigmaDotTot ! internal entropy produciton at time t
      type(rxnNet)  rNet ! stores reactions, free energies, etc calculated in rxnProperties.  Control variables also get set at t and x.
      character (len=20) fmt
      integer i, j, nextRec
      real(dp) cntlTX(ncntl(3))
      real(dp) D, dDdx, a, dadx, q, dqdx, qL, cL(nC), flux(nC)
      real(dp) vS(1:nC)

      if (myRank /= 0) return ! only done by the master process

      ! Get the next record to be written to next.  This will be the same for all direct access files, so just pick one
      inquire(unit=concUnit, nextrec=nextRec)
      ! Increament jTec (number of time points saved) and update tecplot header for all DA files
      jTec = jTec + 1
      write(concUnit    ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(bioSUnit    ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(cntlUnit    ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(rxnUnit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(F_Kunit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(F_Tunit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(dGrUnit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(sigmaRijUnit,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(sigmaUnit   ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(envUnit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec
      write(flxUnit     ,'(a,i4,a,i5,a,i1)',rec=IJKrec) 'I = ', iTec,', J = ', jTec, ', K = ', kTec

      do i=1,nXpts
         associate ( x => xPts(i), conc => Cout(1:nconc,i,1), bioS => Cout(nconc+1:nconc+nbioS,i,1) )   
            ! Interpolate time in tOCnodes and x in xOCnodes, since OC grid is typically coarser than PDE grid
            call interp2D (ncntl(1)+1, ncntl(2), ncntl(3), tOCnodes, xOCnodes, cntl, t, x, cntlTX)
            ! get reaction properties
            call rxnProperties(t, x, conc, bioS, ncntl(3), cntlTX, rNet, sumSigmaDot)
            ! Save outputs at current position x
            fmt = '(   (g15.7,1x))'
            write(fmt(2:4),'(i3.3)') nconc+2
            write(concUnit, fmt, rec=nextRec) t, x, conc(1:nconc)
            
            write(fmt(2:4),'(i3.3)') nbioS+2
            write(bioSUnit, fmt, rec=nextRec) t, x, bioS(1:nbioS)

            ! control variables have been set by rxnProperties now
            write(fmt(2:4),'(i3.3)') 2+34
            write(cntlUnit, fmt, rec=nextRec) t, x, rNet%phyC%eps, rNet%phyC%omg(1:2), rNet%GSBC%eps, rNet%GSBC%omg(1:2), rNet%GzC%eps, rNet%GzC%omg(1:8), &
                           rNet%AGzC%eps, rNet%AGzC%omg(1:8), rNet%BacC%eps, rNet%BacC%omg(1:3), rNet%SRBC%eps, rNet%SRBC%omg(1:3), rNet%PhC%eps, rNet%SOxC%eps
            
            ! Reaction rates, drivers and associated free energies are now stored in rNet, also, two files are now used instead of one
            write(fmt(2:4),'(i3.3)') ncntl(3)+2
            write(rxnUnit , fmt, rec=nextRec) t, x, rNet%PhyR%r(1:2), rNet%GSBR%r(1:2), rNet%GzR%r(1:8), rNet%AGzR%r(1:8), rNet%BacR%r(1:3), rNet%SRBR%r(1:3), &
               rNet%PhR%r(1), rNet%SOxR%r(1)
            ! F_K
            write(F_Kunit , fmt, rec=nextRec) t, x, rNet%PhyR%F_K(1:2), rNet%GSBR%F_K(1:2), rNet%GzR%F_K(1:8), rNet%AGzR%F_K(1:8), rNet%BacR%F_K(1:3), rNet%SRBR%F_K(1:3), &
               rNet%PhR%F_K(1), rNet%SOxR%F_K(1)
            ! F_T
            write(F_Tunit , fmt, rec=nextRec) t, x, rNet%PhyR%F_T(1:2), rNet%GSBR%F_T(1:2), rNet%GzR%F_T(1:8), rNet%AGzR%F_T(1:8), rNet%BacR%F_T(1:3), rNet%SRBR%F_T(1:3), &
               rNet%PhR%F_T(1), rNet%SOxR%F_T(1)
            ! dG
            write(fmt(2:4),'(i3.3)') ncntl(3)+2
            write(dGrUnit , fmt, rec=nextRec) t, x, rNet%PhyR%dG(1:2), rNet%GSBR%dG(1:2), rNet%GzR%dG(1:8), rNet%AGzR%dG(1:8), rNet%BacR%dG(1:3), rNet%SRBR%dG(1:3), &
               rNet%PhR%dG(1), rNet%SOxR%dG(1)

            ! sigmaRij
            write(fmt(2:4),'(i3.3)') ncntl(3)+2
            write(sigmaRijUnit , fmt, rec=nextRec) t, x, rNet%PhyR%sig(1:2), rNet%GSBR%sig(1:2), rNet%GzR%sig(1:8), rNet%AGzR%sig(1:8), rNet%BacR%sig(1:3), rNet%SRBR%sig(1:3), &
               rNet%PhR%sig(1), rNet%SOxR%sig(1)
            
            sumSigmaDotTot = sum(sumSigmaDot)
            write(sigmaUnit, '(7(g15.7,1x))', rec=nextRec) t, x, sumSigmaDotTot, sumSigmaX(i), &
               100._dp*sumSigmaDot(1)/sumSigmaDotTot, 100._dp* sumSigmaDot(2)/sumSigmaDotTot,  100._dp*sumSigmaDot(3)/sumSigmaDotTot
            
            write(fmt(2:4),'(i3.3)') 5
            write(envUnit,fmt, rec=nextRec) t, x, I0(t)*exp( kIntg_x (x) ), I0(t)*exp( kIntg_x (x) )*1000._dp/86400._dp, kIntg_x (x)
            
            call dispersion (t, x, D, dDdx) ! Dispersion and it's derivative (m2/d)
            call ccArea     (t, x, a, dadx) ! crossectional area (m2)
            call volFlow    (t, x, q, dqdx, qL, cL) ! volumetric flow rate (m3/d) and laterial inputs if any. qL is in m2/d
            ! calclate flux of constituent assocaite with advection and dispersion, full dimensions
            vS = vSink*stepUp(t, tStep, tSig) ! this slowly spins up sinking.
            flux(1:nC) = - D*a*Cout(1:nC,i,2) + q*Cout(1:nC,i,1) + vS(1:nC)*a*Cout(1:nC,i,1) 
            write(fmt(2:4),'(i3.3)') nC+2
            write(flxUnit, fmt, rec=nextRec) t, x, flux(1:nC)            
         end associate
         nextRec = nextRec + 1
      end do
      ! Save BACOLI info
      write(BACunit,'(g15.7,6(1x,i7))') t, ADRsol%nint, ADRsol%num_remeshings, ADRsol%num_ini_remeshings, ADRsol%num_cold_restarts, ADRsol%num_accepted_time_steps, ADRsol%prev_bdf_order
      return
   end subroutine saveSolution
   
   subroutine saveGrid (ADRsol)
      ! This saves the grid points
      use bacoli95_mod, only: bacoli95_sol 
      use parameters, only: Xc, gridUnit, myRank
      implicit none
      integer i
      type(bacoli95_sol), intent(in):: ADRsol 
      if (myRank /= 0) return ! only master process to save output
      ! BACOLI grid
      do i=1,ADRsol%nint+1
         write(gridUnit,'(2(g15.7,1x))') ADRsol%t0, ADRsol%x(i)*Xc
      end do
      return
   end subroutine saveGrid
   
   subroutine saveOCgrid (tIntf)      
      ! this routine save the xOC and tOC grid points over an interval
      use realPrec
      use parameters      
      implicit none
      real(dp), intent(in):: tIntf ! time (d) of current end of interval.
      ! local declarations
      integer i, j
      if (myRank /= 0) return ! only the master rank saves data.
      do i=0,ncntl(1)
         do j=1,ncntl(2)
            if (tOCnodes(i) > tIntf) return ! only save nodes within the finite interval
            write(OCgridUnit,'(2(g15.7,1x))') tOCnodes(i), xOCnodes(j)
         end do
      end do
      return
   end subroutine saveOCgrid   

   subroutine summarySave(i, tInt0, sumSigmaTotInf, sumSigmaTotFin)
      ! This routine saves the interval information and displays to screen
      use realPrec
      use parameters
      use mpi, only: MPI_Wtime
      implicit none
      integer, intent(in):: i ! current interval
      real(dp), intent(in):: tInt0 ! End time of current interval
      real(dp), intent(in):: sumSigmaTotInf ! total entropy production over infinate interval integrated over Xpts (J/K)
      real(dp), intent(in):: sumSigmaTotFin ! total entropy production over finite interval integrated over Xpts (J/K)
      ! local declarations
      real(dp) wallTime
      
      if (myRank /= 0) return ! only the master process is allowed to save data
      ! get time
      WTime_f = MPI_Wtime() ! end timer in seconds
      wallTime = (WTime_f-WTime_0)/60._dp
      write(summaryUnit,'(i4,1x,2(g15.7,1x),i5,3(g15.7,1x),i7)') i, tInt0, wallTime, localIter, localRes, sumSigmaTotInf, sumSigmaTotFin, errCntTot 
      ! Display output to screen
      write(*,'(a)') 'Interval  tIntf  CPU(min)   sigmaInf        sigmaInt      errCnt'                                                                                         
      write(*,'(1x,i4,5x,f5.1,2x,f8.1,1x,g13.5,3x,g13.5,1x,i6/)') i, tInt0,  wallTime, sumSigmaTotInf, sumSigmaTotFin, errCntTot
      return
   end subroutine summarySave   
   
   subroutine copySol (solOrig, solCopy)
      ! this routine is used to copy the VALUES, not pointers of solOrg to solCopy.
      ! if you just use solCopy = solOrig,  only pointers are copied.
      use bacoli95_mod, only: bacoli95_sol
      implicit none
      type(bacoli95_sol):: solOrig, solCopy
      
      ! begin the copy
      solCopy%npde = solOrig%npde
      solCopy%nint_max = solOrig%nint_max
      solCopy%kcol = solOrig%kcol
      solCopy%estimator = solOrig%estimator
      solCopy%maxord = solOrig%maxord
      solCopy%atol = solOrig%atol
      solCopy%rtol = solOrig%rtol
      solCopy%ini_ss = solOrig%ini_ss
      solCopy%t0 = solOrig%t0
      solCopy%x = solOrig%x
      solCopy%nint = solOrig%nint
      solCopy%mflag = solOrig%mflag
      solCopy%idid = solOrig%idid
      solCopy%y = solOrig%y
      solCopy%rpar = solOrig%rpar
      solCopy%ipar = solOrig%ipar
      solCopy%num_remeshings = solOrig%num_remeshings
      solCopy%num_ini_remeshings = solOrig%num_ini_remeshings
      solCopy%num_cold_restarts = solOrig%num_cold_restarts
      solCopy%num_accepted_time_steps = solOrig%num_accepted_time_steps
      solCopy%min_len_ipar = solOrig%min_len_ipar
      solCopy%min_len_rpar = solOrig%min_len_rpar
      solCopy%prev_bdf_order = solOrig%prev_bdf_order
      solCopy%prev_time_step_size = solOrig%prev_time_step_size
      return
   end subroutine copySol
   
   subroutine getADRsol (iostat)
      ! This gets the ADRsol from file if readADRsol is .true.
      use realPrec
      use parameters
      use bacoli95_mod, only: bacoli95_sol
      integer, intent(out)           :: iostat ! if an error occurs reading ADRsol
      ! local declarations
      integer i, unitNo
      
      solGood = .false.
      iostat = 0
      if (.not. readADRsol) return
      ! open the t0 sol file.  Make this read only for MPI
      open(newunit=unitNo, file=trim(basefilename)//'_t0.sol', form='unformatted', status='old', action='read')
      call read_sol(unitNo, ADRsol, iostat)
      close(unit=unitNo)
      if (iostat /= 0) return
      solGood = .true.      
      return
   end subroutine getADRsol   
   
   subroutine write_sol(unit, sol, iostat)
      ! used to write sol of BACOLI95 to a file.
      ! unit must alreadyt be open, such as:
      ! open(newunit=solUnit, file='test.sol', form='unformatted', status='replace', action='write')   
      use bacoli95_mod, only: bacoli95_sol
      implicit none
      integer, intent(in)             :: unit
      type(bacoli95_sol), intent(in)  :: sol
      integer, intent(out)            :: iostat

      ! Write a record for each var, but add sizes to those that are allocatable
      write(unit, iostat=iostat) sol%npde
      write(unit, iostat=iostat) sol%nint_max
      write(unit, iostat=iostat) sol%kcol
      write(unit, iostat=iostat) sol%estimator
      write(unit, iostat=iostat) sol%maxord
      write(unit, iostat=iostat) size(sol%atol)
      write(unit, iostat=iostat) sol%atol
      write(unit, iostat=iostat) size(sol%rtol)
      write(unit, iostat=iostat) sol%rtol
      write(unit, iostat=iostat) sol%ini_ss
      write(unit, iostat=iostat) sol%t0
      write(unit, iostat=iostat) size(sol%x)
      write(unit, iostat=iostat) sol%x
      write(unit, iostat=iostat) sol%nint
      write(unit, iostat=iostat) sol%mflag
      write(unit, iostat=iostat) sol%idid
      write(unit, iostat=iostat) size(sol%y)
      write(unit, iostat=iostat) sol%y
      write(unit, iostat=iostat) size(sol%rpar)
      write(unit, iostat=iostat) sol%rpar
      write(unit, iostat=iostat) size(sol%ipar)
      write(unit, iostat=iostat) sol%ipar
      write(unit, iostat=iostat) sol%num_remeshings
      write(unit, iostat=iostat) sol%num_ini_remeshings
      write(unit, iostat=iostat) sol%num_cold_restarts
      write(unit, iostat=iostat) sol%num_accepted_time_steps
      write(unit, iostat=iostat) sol%min_len_ipar
      write(unit, iostat=iostat) sol%min_len_rpar
      write(unit, iostat=iostat) sol%prev_bdf_order
      write(unit, iostat=iostat) sol%prev_time_step_size  
      return
   end subroutine write_sol

   subroutine read_sol(unit, sol, iostat)
      ! used to read a stored sol of BACOLI95
      ! unit should be opened like: open(newunit=solUnit, file=trim(filename)//'.sol', form='unformatted', status='old', action='read')
      ! although, specifying action='read' is probably overkill
      use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown
      implicit none
      integer, intent(in)             :: unit  ! Unit connected to file containing sol.
      type(bacoli95_sol), intent(out) :: sol
      integer, intent(out)            :: iostat

      integer allocSize
      ! read is sol type, but allocate pointer space as needed, as it is ASSUMED not to be allocated
      read(unit, iostat=iostat) sol%npde
      read(unit, iostat=iostat) sol%nint_max
      read(unit, iostat=iostat) sol%kcol
      read(unit, iostat=iostat) sol%estimator
      read(unit, iostat=iostat) sol%maxord
      read(unit, iostat=iostat) allocSize
      allocate( sol%atol(allocSize) )      
      read(unit, iostat=iostat) sol%atol
      read(unit, iostat=iostat) allocSize
      allocate( sol%rtol(allocSize) ) 
      read(unit, iostat=iostat) sol%rtol
      read(unit, iostat=iostat) sol%ini_ss
      read(unit, iostat=iostat) sol%t0
      read(unit, iostat=iostat) allocSize
      allocate( sol%x(allocSize) ) 
      read(unit, iostat=iostat) sol%x
      read(unit, iostat=iostat) sol%nint
      read(unit, iostat=iostat) sol%mflag
      read(unit, iostat=iostat) sol%idid
      read(unit, iostat=iostat) allocSize
      allocate( sol%y(allocSize) ) 
      read(unit, iostat=iostat) sol%y
      read(unit, iostat=iostat) allocSize
      allocate( sol%rpar(allocSize) ) 
      read(unit, iostat=iostat) sol%rpar
      read(unit, iostat=iostat) allocSize
      allocate( sol%ipar(allocSize) ) 
      read(unit, iostat=iostat) sol%ipar
      read(unit, iostat=iostat) sol%num_remeshings
      read(unit, iostat=iostat) sol%num_ini_remeshings
      read(unit, iostat=iostat) sol%num_cold_restarts
      read(unit, iostat=iostat) sol%num_accepted_time_steps
      read(unit, iostat=iostat) sol%min_len_ipar
      read(unit, iostat=iostat) sol%min_len_rpar
      read(unit, iostat=iostat) sol%prev_bdf_order
      read(unit, iostat=iostat) sol%prev_time_step_size   
      return
   end subroutine read_sol             
   
   subroutine setCnames(conc, bioS, C)
      ! This routine passes conc and bioS values to useful names
      use realPrec
      use typeVars, only: concentrations
      use parameters, only: absZero, nconc, nbioS, nh3_hold, N_D_hold
      use functions, only: fZero ! This is used to prevent state variable from being <= 0.  This transition to absZero occurs at 2*absZero
      implicit none
      real(dp),             intent(in):: conc(nconc) ! vector of concentration variables
      real(dp),             intent(in):: bioS(nbioS) ! vector of biological structure variables
      type(concentrations), intent(out):: C          ! derived type that hold state names. See also module typeVars
      ! Begin assignments
      C%sal   = fZero(absZero, conc(1 ))
      C%o2    = fZero(absZero, conc(2 ))
      C%dic   = fZero(absZero, conc(3 ))
      C%nh3   = nh3_hold ! NH3 is not currently a state var, but it is used in delG calculations, etc
      C%po4   = fZero(absZero, conc(4 ))
      C%so4   = fZero(absZero, conc(5 ))
      C%h2s   = fZero(absZero, conc(6 ))
      C%C_Phy = fZero(absZero, conc(7 ))
      C%C_GSB = fZero(absZero, conc(8 ))
      C%C_L   = fZero(absZero, conc(9 ))
      C%C_D   = fZero(absZero, conc(10))
      C%N_D   = N_D_hold ! N_D (detrital N) is not currently a state var, but it is used in delG calculations, etc
      C%P_D   = fZero(absZero, conc(11))
      C%PhyS  = fZero(absZero, bioS(1 ))
      C%GSBS  = fZero(absZero, bioS(2 ))
      C%GzS   = fZero(absZero, bioS(3 ))
      C%AGzS  = fZero(absZero, bioS(4 ))
      C%BacS  = fZero(absZero, bioS(5 ))
      C%SRBS  = fZero(absZero, bioS(6 ))
      C%PhS   = fZero(absZero, bioS(7 ))
      C%SOxS  = fZero(absZero, bioS(8 ))
      return
   end subroutine setCnames   
   
   subroutine rxnProperties(t, x, conc, bioS, nU, cntlTX, rNet, sumSigmaDot)
      ! Calculates reaction rates, and Gibbs free energy of rxn, and total entropy production by reactions
      ! THIS ROUTINE WILL DEPEND OF THE PROBLEM BEING SOLVED
      ! This takes a modular-like call, in that each metabolic sub group has it's own subroutine,
      ! which should make adding or removing functional groups slightly easier.
      ! For now, do not implment omega partitioning to sub reactions in group
      ! This routine can be called from outside of the ODE routine that mainly uses it.  This is useful
      ! for saving rxn and free energy over an interval during storage.
      ! Reactions are evaluated at time t and at location x only.  
      ! The state variables are (in order saved in conc and bioS)
      ! conc: sal, o2, dic, NH3, H3PO4, C_p, C_D, N_D, P_D  (9 currently)
      use realPrec
      use parameters
      use typeVars
      use thermoData
      use functions, only: I0      
      use interfaces
      implicit none
      real(dp),     intent(in) :: t                ! current time (d)
      real(dp),     intent(in) :: x                ! current location (m)
      real(dp),     intent(in) :: conc(nconc)      ! chemical concentration state variables at (t,x)
      real(dp),     intent(in) :: bioS(nbioS)      ! concentration of biological structures at (t,x)
      integer,      intent(in) :: nU               ! number of control variables
      real(dp),     intent(in) :: cntlTX(nU)       ! This assumes that only the control variables at the current time and location are passed.
      type(rxnNet), intent(out):: rNet             ! This contains all the reaction coef as well as control variables that will be set below
      real(dp),     intent(out):: sumSigmaDot(3)   ! entropy production summed over all reactions (1), water (2) and particles (3), per unit length (J/m/d/K)
      ! Local declarations
      real(dp) dadx
      type(concentrations) C ! used to pass concentrations to subroutines
      type (rxnVariables):: rV ! see module rxnVars for list of variables and parameters that need to be set here. 
      real(dp) sumSigmaDotRxn, sumSigmaDotH2O, sumSigmaDotPart
      
      ! calculations used by rxn routines
      ! set names, and insure none are <= to absZero
      call setCnames (conc, bioS, C)
      
      ! setup the reaction variables derived type (rV here)
      call ccArea (t, x, rV%csAtx, dadx) ! get cross sectional area, deravative not used here.
      call tempK  (t, x, rV%T_K) ! get temperature in K
      call pHtx   (t, x, rV%pH ) ! get pH
      rV%is          = 0.72*C%sal/35.0 ! use the ionic strength based on a seawater value of 0.7 M at (see http://www.teos-10.org/pubs/gsw/html/gsw_ionic_strength_from_SA.html).
      rV%absZero     = absZero ! as defined in parameters module
      rV%depth       = depth
      ! calculate necessary thermodynamic variables, such as standard free energy's of formation (since many routines will need)  See thermoData module
      rV%dGf0_bioS   = dGf0(yeastsp,   rV%T_K,rV%is,rV%pH)  ! free energies of formation of species at current temp (K), is (M), and pH. (kJ/mol)
      rV%dGf0_C_D    = dGf0(ch2osp,    rV%T_K,rV%is,rV%pH)  ! Detritus C (use ch2o), N (use nh3) and P (use h3po4).
      rV%dGf0_N_D    = dGf0(ammoniasp, rV%T_K,rV%is,rV%pH)
      rV%dGf0_P_D    = dGf0(pisp,      rV%T_K,rV%is,rV%pH)
      rV%dGf0_ch2o   = dGf0(ch2osp,    rV%T_K,rV%is,rV%pH)   ! glucose, single carbon
      !rV%dGf0_h2aq   = dGf0(h2aqsp,    rV%T_K,rV%is,rV%pH)  ! Not using yet.
      rV%dGf0_o2aq   = dGf0(o2aqsp,    rV%T_K,rV%is,rV%pH)
      rV%dGf0_h2so4  = dGf0(sulfatesp, rV%T_K,rV%is,rV%pH)
      rV%dGf0_h2saq  = dGf0(h2saqsp,   rV%T_K,rV%is,rV%pH)
      rV%dGf0_dic    = dGf0(co2totsp,  rV%T_K,rV%is,rV%pH) ! or DIC
      !rV%dGf0_hno3   = dGf0(nitratesp, rV%T_K,rV%is,rV%pH) ! Not using yet
      rV%dGf0_nh3    = dGf0(ammoniasp, rV%T_K,rV%is,rV%pH)
      rV%dGf0_h3po4  = dGf0(pisp,      rV%T_K,rV%is,rV%pH)
      rV%dGf0_h2o    = dGf0(h2osp,rV%T_K,rV%is,rV%pH)
      rV%dGr_Ggamma  = delrGgamma  ! Energy available in green light, accounting for thermo efficiency (J/mmol) set in input file
      rV%co2         = freeCO2 (rV%T_K,rV%is,rV%pH)*C%dic ! Free CO2 concentration in mmol/m3      
      rV%hco3        = freeHCO3(rV%T_K,rV%is,rV%pH)*C%dic ! Free HCO3 concentration in mmol/m3      
      rV%RkJ         = RkJ         ! gas constant in J/K/mmol or kJ/K/mol (this parameter set in thermoData module)
      rV%nuStar      = nuStar      ! max specific growth rate (1/d)
      rV%nuDet       = nuDet       ! max specific growth rate on detritus (1/d)
      rV%kappa       = kappa       ! substrate affinity (mmol/m3)
      rV%alp         = bioS_H      ! H in biological structures (all the same for now, and single C basis)
      rV%bet         = bioS_O      ! O in biological structures (all the same for now, and single C basis)
      rV%gam         = bioS_N      ! N in biological structures (all the same for now, and single C basis)
      rV%del         = bioS_P      ! P in biological structures (all the same for now, and single C basis)
      rV%k_w         = k_w         ! light absorption by water (1/m)
      rV%k_p         = k_p         ! light absoprtion coef by algae and other particles (m^2/mmol-S)
      rV%cell_F      = cell_F      ! The true concentration of C_P is much higher that the state solution for because it is only inside the cell.  Adjust for kinetics
      rV%Itx         = I0(t)*exp( kIntg_x (x) ) ! light level at (t,x) (mmol photons /m^2 /d (Not micromoles)) 
   
      sumSigmaDotRxn = 0_dp
      ! Get reactions associated with all the Biological Structures. 
      ! Phy Biological Structure
      rNet%phyC%eps = cntlTX(1) ! eps for Phy
      call omega_free (2, cntlTX(2:2), rNet%phyC%omg) ! converts one w's to two Omega's
      call Phy_bioS (t, x, C, rNet%phyC, rV, rNet%phyR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%phyR%sig(1:2)) ! add entropy to total
      
      ! GSB Biological Structure
      rNet%GSBC%eps = cntlTX(3) ! eps for GSB
      call omega_free (2, cntlTX(4:4), rNet%GSBC%omg) ! converts one w's to two Omega's
      call GSB_bioS (t, x, C, rNet%GSBC, rV, rNet%GSBR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%GSBR%sig(1:2)) ! add entropy to total
      
      ! Gz Biological Structure
      rNet%GzC%eps = cntlTX(5) ! eps for Gz
      call omega_free (8, cntlTX(6:12), rNet%GzC%omg) ! converts w's to Omega's
      call Gz_bioS (t, x, C, rNet%GzC, rV, rNet%GzR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%GzR%sig(1:8)) ! add entropy to total
      
      ! AGz Biological Structure
      rNet%AGzC%eps = cntlTX(13) ! eps for AGz
      call omega_free (8, cntlTX(14:20), rNet%AGzC%omg) ! converts w's to Omega's
      call AGz_bioS (t, x, C, rNet%AGzC, rV, rNet%AGzR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%AGzR%sig(1:8)) ! add entropy to total
 
      ! Bac Biological Structure
      rNet%BacC%eps = cntlTX(21) ! eps for Bac
      call omega_free (3, cntlTX(22:23), rNet%BacC%omg) ! converts w's to Omega's
      call Bac_bioS (t, x, C, rNet%BacC, rV, rNet%BacR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%BacR%sig(1:3)) ! add entropy to total

      ! SRB Biological Structure
      rNet%SRBC%eps = cntlTX(24) ! eps for SRB
      call omega_free (3, cntlTX(25:26), rNet%SRBC%omg) ! converts w's to Omega's
      call SRB_bioS (t, x, C, rNet%SRBC, rV, rNet%SRBR)
      sumSigmaDotRxn = sumSigmaDotRxn + sum(rNet%SRBR%sig(1:3)) ! add entropy to total
      
      ! Ph Biological Structure
      rNet%PhC%eps = cntlTX(27) ! eps for SRB
      ! Only omg for Ph is 1 by definition
      call Ph_bioS (t, x, C, rNet%PhC, rV, rNet%PhR)
      sumSigmaDotRxn = sumSigmaDotRxn + rNet%PhR%sig(1) ! add entropy to total
      
      ! SOx Sulfur bacteria
      rNet%SOxC%eps = cntlTX(28) ! There is only one control variable for SOx
      call SOx_bioS (t, x, C, rNet%SOxC, rV, rNet%SOxR)
      sumSigmaDotRxn = sumSigmaDotRxn + rNet%SOxR%sig(1) ! add entropy to total

      ! Add the entropy associated with light dissipation on all other particles, water, etc. and contribuion from the non photoactive components
      ! of phototrophs (the photo active compenents have already be accounted for).
      ! It is also assume that carbon stored in phototrophs (C_Phy and C_GSB) also intercepts and dissipates light.
      ! Note, entropy associated with light adsorption for photoheterotrophs, Ph, has already be included above. 
      ! Entropy production associated with water absorption
      sumSigmaDotH2O = - rV%Itx*rV%dGr_Ggamma*rV%k_w*rV%csAtx/rV%T_K
      ! Entropy production associated with non-photosynthetic parts of particles
      sumSigmaDotPart = - rV%Itx*rV%dGr_Ggamma*rV%k_p*( rNet%phyC%omg(2)*C%PhyS + C%C_Phy + rNet%GSBC%omg(2)*C%GSBS + C%C_GSB &
                  & + C%C_D + C%GzS + C%AGzS + C%BacS + C%SRBS + C%SOxS )*rV%csAtx/rV%T_K
      ! Total entropy production
      ! sumSigmaDotTot = sumSigmaDotRxn + sumSigmaDotH2O + sumSigmaDotPart
      sumSigmaDot(1) = sumSigmaDotRxn
      sumSigmaDot(2) = sumSigmaDotH2O
      sumSigmaDot(3) = sumSigmaDotPart
      return
   end subroutine rxnProperties
   
   subroutine Phy_bioS (t, x, C, phyC, rV, phyR)
      ! This routine is use to model aerobic phytoplankton (1: Phy)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! An alternative would be in add a derived type that could package all the parameters, but that would be for later development
      ! This version does not include a mortality reaction, so only co2 fixation to ch2o (cL) and converstion of that to phyS
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a) and biological Structure conc (conc. in mmol/m3)
      type(PhyCntl)       , intent(in) :: phyC  ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(PhyRxns)       , intent(out):: phyR  ! all outputs associated with the aerobic phytoplankton reactions.  See module typeVars, type PhyRxns
      
      ! local declarations
      real(dp) delI_P ! Light captured for photosynthesis by the component of phytoplankton allocated to photosynthesis, given by omg(1)
      real(dp) dG0_r1, dG_r1 ! free energy for reaction 1, CO2 fixation (kJ/mol or J/mmol)
      real(dp) dG0_r2A, dG_r2A, dG0_r2C, dG_r2C
      real(dp) concFac ! testing a concentration factor
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Begin calculations
      ! Reaction 1, CO2 fixation
      ! eps h2co3 + n(1) hv -> eps(C_P + o2)  NOTE, now considering H2CO3 to include CO2 and biocarbonate.
      ! Calculate light captured by Phy photosynthetic apparatus only. 
      delI_P  = rV%k_p*phyC%omg(1)*C%phyS*rV%Itx ! note, only omg(1) fraction of phyS is allocated to light harvesting for chemosynthesis.
      ! calculate free energy associated with CH2O3 -> CH2O + O2
      dG0_r1 = (rV%dGf0_ch2o + rV%dGf0_o2aq) - rV%dGf0_dic ! CO2 + H2O fixation into glucose and oxygen (unit carbon)   
      dG_r1  = dG0_r1 + rV%RkJ*rV%T_K*( (log(C%C_Phy)+ln106) + (log(C%o2)+ln106) - (log(C%dic)+ln106)  )
      ! calculate the mmol of photons needed to fix 1 mmol co2
      phyR%n(1) = -dG_r1/rV%dGr_Ggamma ! this is the n1 coef.
      phyR%dG(1) = -(1_dp - phyC%eps)*dG_r1 ! free energy of reaction for reaction 1, CO2 fixation
      phyR%ne(1) = phyR%n(1) ! electrons transfered in catabolic reaction (in version 1.1, this was fixed at 4)
      PhyR%F_T(1) = F_Thermo(phyR%dG(1), phyR%ne(1), rV%T_K)
      PhyR%F_K(1) =  (rV%co2+rV%hco3)/(rV%co2+rV%hco3 + rV%kappa*phyC%eps**4) ! added bicarbonate uptake to kinetics
      phyR%r(1) = delI_P/phyR%n(1)*PhyR%F_T(1)*PhyR%F_K(1)
      phyR%sig(1) = dG_r1*( delI_P/phyR%n(1) - phyC%eps*phyR%r(1) )*rV%csAtx/rV%T_K  ! EP from light dissipation by reaction 1 (per unit length for 1D problem)
         
      ! Calculate for phototroph biosynthesis reaction, r2
      ! (1+eps n2)C_Phy + eps(gam NH3 + del H2PO4) + [1 + eps(aA2 + n2 -1)]O2 -> eps PhyS + eps bA2 H2O + (1 + eps(n2-1))dic
      ! calculate reaction coefficients, but only anabolic reactions have them
      phyR%aA2 = (      - rV%alp + 2._dp*rV%bet + 3._dp*rV%gam - 5._dp*rV%del)/4._dp
      phyR%bA2 = (2._dp - rV%alp                + 3._dp*rV%gam + 3._dp*rV%del)/2._dp
      ! anabolic reaction thermodynamics
      dG0_r2A = (rV%dGf0_bioS +  phyR%bA2*rV%dGf0_h2o) - (rV%dGf0_ch2o + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4 + phyR%aA2*rV%dGf0_o2aq)
      dG_r2A  = dG0_r2A + rV%RkJ*rV%T_K*( (log(C%phyS)+ln106) - (log(C%C_Phy)+ln106) - rV%gam*(log(C%nh3)+ln106) - rV%del*(log(C%po4)+ln106) - phyR%aA2*(log(C%o2)+ln106)  )
      ! catabolic reaction thermodynamics
      dG0_r2C = rV%dGf0_dic - ( rV%dGf0_ch2o + rV%dGf0_o2aq )
      dG_r2C = dG0_r2C + rV%RkJ*rV%T_K*( (log(C%dic)+ln106) - (log(C%C_Phy)+ln106) - (log(C%o2)+ln106) )
      ! coupling of anabolic to catabolic reaction
      phyR%n(2) = - dG_r2A/dG_r2C
      ! Free energy of combined whole reaction
      phyR%dG(2) = (1._dp - phyC%eps)*dG_r2C
      ! Reaction (2) rate
      phyR%ne(2) = 4. ! electrons transfered in catabolic reaction
      PhyR%F_T(2) = F_Thermo(phyR%dG(2), phyR%ne(2), rV%T_K)
      PhyR%F_K(2) = ( C%C_Phy*rV%cell_F/(C%C_Phy*rV%cell_F + rV%kappa*phyC%eps**4) )*( C%po4/rV%del/(C%po4/rV%del + rV%kappa*phyC%eps**4)  )*( C%o2/(C%o2 + rV%kappa*phyC%eps**4)  ) ! don't include NH3 limatations yet
      phyR%r(2) = rV%nuStar*phyC%omg(2)*phyC%eps**2*C%phyS*PhyR%F_T(2)*PhyR%F_K(2)
      ! Entropy production
      phyR%sig(2) = - phyR%r(2)*phyR%dG(2)*rV%csAtx/rV%T_K
      return
   end subroutine Phy_bioS
   
   subroutine GSB_bioS (t, x, C, GSBC, rV, GSBR)
      ! This routine is use to model anaerobic photoautotrophs, green sulfur bacteria (2: GSB)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! An alternative would be in add a derived type that could package all the parameters, but that would be for later development
      ! Two reactions are considered here: co2 fixation to ch2o (cL) and converstion of that to GSBS
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a) and biological Structure conc (conc. in mmol/m3)
      type(GSBCntl)       , intent(in) :: GSBC  ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(GSBRxns)       , intent(out):: GSBR  ! all outputs associated with the GSB reactions.  See module typeVars, type GSBRxns
      
      ! local declarations
      real(dp) delI_GSB ! Light captured for photosynthesis by the component of GSB allocated to photosynthesis, given by omg(1)
      real(dp) dG0_r1, dG_r1 ! free energy for reaction 1, CO2 fixation (kJ/mol or J/mmol)
      real(dp) dG0_r2A, dG_r2A, dG0_r2C, dG_r2C
      real(dp) concFac ! testing a concentration factor
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Begin calculations
      ! Reaction 1, CO2 fixation
      ! eps (h2co3 + 1/2H2S) + n(1) hv -> eps(C_GSB + 1/2H2SO4)  NOTE, now considering H2CO3 to include CO2 and biocarbonate.
      ! Calculate light captured by GSB photosynthetic apparatus only. 
      delI_GSB  = rV%k_p*GSBC%omg(1)*C%GSBS*rV%Itx ! note, only omg(1) fraction of GSBS is allocated to light harvesting for chemosynthesis.
      ! calculate free energy associated with CH2O3 -> CH2O + O2
      dG0_r1 = (rV%dGf0_ch2o + rV%dGf0_h2so4/2._dp) - (rV%dGf0_dic + rV%dGf0_h2saq/2._dp) ! CO2 + H2S fixation into glucose and H2SO4 (unit carbon)   
      dG_r1  = dG0_r1 + rV%RkJ*rV%T_K*( (log(C%C_GSB)+ln106) + (log(C%so4)+ln106)/2._dp - (log(C%dic)+ln106) - (log(C%h2s)+ln106)/2._dp )
      ! calculate the mmol of photons needed to fix 1 mmol co2
      GSBR%n(1) = max(-dG_r1/rV%dGr_Ggamma, 0.1_dp) ! this is the n1 coef., set a minimum value of 1 since dG_r1 can sometimes be <0
      GSBR%dG(1) = GSBC%eps*dG_r1 + GSBR%n(1)*rV%dGr_Ggamma !-(1_dp - GSBC%eps)*dG_r1 ! free energy of reaction for reaction 1, CO2 fixation
      GSBR%ne(1) = GSBR%n(1) ! electrons transfered in catabolic reaction.  In V1.1, this was set to 4.  See pg. 83 of notes.
      GSBR%F_T(1) = F_Thermo(GSBR%dG(1), GSBR%ne(1), rV%T_K)
      GSBR%F_K(1) =  (rV%co2+rV%hco3)/(rV%co2+rV%hco3 + rV%kappa*GSBC%eps**4)*( 2._dp*C%h2s/(2._dp*C%h2s + rV%kappa*GSBC%eps**4) ) ! added bicarbonate uptake to kinetics
      GSBR%r(1) = delI_GSB/GSBR%n(1)*GSBR%F_K(1) * GSBR%F_T(1)
      GSBR%sig(1) = -( delI_GSB*rV%dGr_Ggamma + GSBC%eps*dG_r1*GSBR%r(1) )*rV%csAtx/rV%T_K  ! EP from light dissipation by reaction 1 (per unit length for 1D problem)
         
      ! Calculate for phototroph biosynthesis reaction, r2
      ! (1+eps n2)C_GSB + eps(gam NH3 + del H2PO4) + [1/2 + eps(aA2 + n2/2 -1/2)]H2SO4 -> eps GSBP + eps bA2 H2O + [1 + eps(n2-1)]dic + [1/2 + eps(aA2 + n2/2 -1/2)]H2S
      ! calculate reaction coefficients, but only anabolic reactions have them
      GSBR%aA2 = (      - rV%alp + 2._dp*rV%bet + 3._dp*rV%gam - 5._dp*rV%del)/8._dp
      GSBR%bA2 = (2._dp - rV%alp                + 3._dp*rV%gam + 3._dp*rV%del)/2._dp
      ! anabolic reaction thermodynamics
      dG0_r2A = (rV%dGf0_bioS +  GSBR%bA2*rV%dGf0_h2o + GSBR%aA2*rV%dGf0_h2saq) - (rV%dGf0_ch2o + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4 + GSBR%aA2*rV%dGf0_h2so4)
      dG_r2A  = dG0_r2A + rV%RkJ*rV%T_K*( (log(C%GSBS)+ln106) + GSBR%aA2*(log(C%h2s)+ln106) - (log(C%C_GSB)+ln106) - rV%gam*(log(C%nh3)+ln106) &
              & - rV%del*(log(C%po4)+ln106) - GSBR%aA2*(log(C%so4)+ln106)  )
      ! catabolic reaction thermodynamics
      dG0_r2C = rV%dGf0_dic + rV%dGf0_h2saq/2._dp - ( rV%dGf0_ch2o + rV%dGf0_h2so4/2._dp )
      dG_r2C = dG0_r2C + rV%RkJ*rV%T_K*( (log(C%dic)+ln106) + (log(C%h2s)+ln106)/2._dp - (log(C%C_GSB)+ln106) - (log(C%so4)+ln106)/2._dp )
      ! coupling of anabolic to catabolic reaction
      GSBR%n(2) = - dG_r2A/dG_r2C
      ! Free energy of combined whole reaction
      GSBR%dG(2) = (1._dp - GSBC%eps)*dG_r2C
      ! Reaction (2) rate
      GSBR%ne(2) = 4. ! electrons transfered in catabolic reaction
      GSBR%F_T(2) = F_Thermo(GSBR%dG(2), GSBR%ne(2), rV%T_K)
      GSBR%F_K(2) = ( C%C_GSB*rV%cell_F/(C%C_GSB*rV%cell_F + rV%kappa*GSBC%eps**4) )*( C%po4/rV%del/(C%po4/rV%del + rV%kappa*GSBC%eps**4)  )*( C%so4/(C%so4 + rV%kappa*GSBC%eps**4)  ) ! don't include NH3 limatations yet
      GSBR%r(2) = rV%nuStar*GSBC%eps**2*GSBC%omg(2)*C%GSBS*GSBR%F_T(2)*GSBR%F_K(2)
      ! Entropy production
      GSBR%sig(2) = - GSBR%r(2)*GSBR%dG(2)*rV%csAtx/rV%T_K
      return
   end subroutine GSB_bioS
   
   subroutine Gz_bioS (t, x, C, GzC, rV, GzR)
      ! This routine is use to model aerobic grazers (3: Gz)
      ! See Notes for Frontiers MS
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a, etc) and biological Structure conc (conc. in mmol/m3)
      type(GzCntl)        , intent(in) :: GzC   ! Control variables, eps, omg(1:7) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(GzRxns)        , intent(out):: GzR   ! all outputs associated with the aerobic grazer reactions.  See module typeVars
      
      ! local declarations
      integer i
      real(dp) dG0_ri, dG_ri ! Free energy for reaction i, which is the consumption rate of BioS(i) (kJ/mol or J/mmol)
      real(dp) ri ! for most of the reaction rates
      real(dp) F_o2
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
      
      ! Catabolic stiochiometric coef.
      GzR%aCi = (  4._dp + rV%alp - 2._dp*rV%bet - 3._dp*rV%gam + 5._dp*rV%del - 4._dp*GzC%eps)/4._dp
      GzR%bCi = (- 2._dp + rV%alp                - 3._dp*rV%gam - 3._dp*rV%del                )/2._dp
      
      ! Free energy of reactions (Note, this does not depend on C_Phy or C_GSB), so these are all the same
      dG0_ri = GzC%eps*rV%dGf0_bioS + (1._dp - GzC%eps)**2*rV%dGf0_dic + GzC%eps*(1._dp - GzC%eps)*rV%dGf0_ch2o + rv%del*(1._dp - GzC%eps)*rV%dGf0_h3po4 &
             & + rv%gam*(1._dp - GzC%eps)*rV%dGf0_nh3 + GzR%bCi*(1._dp - GzC%eps)*rV%dGf0_h2o - ( rV%dGf0_bioS + GzR%aCi*(1._dp - GzC%eps)*rV%dGf0_o2aq )
      ! First calculate the quantity that is the same for all reactions.
      dG_ri  = dG0_ri + rV%RkJ*rV%T_K*( GzC%eps*(log(C%GzS)+ln106) + (1._dp - GzC%eps)**2*(log(C%dic)+ln106) + GzC%eps*(1._dp - GzC%eps)*(log(C%C_D)+ln106) &
             + rv%del*(1._dp - GzC%eps)*((1._dp - GzC%eps)*(log(C%po4)+ln106) + GzC%eps*(log(C%P_D)+ln106)) &
             + rv%gam*(1._dp - GzC%eps)*((1._dp - GzC%eps)*(log(C%nh3)+ln106) + GzC%eps*(log(C%N_D)+ln106)) &
             - (GzR%aCi*(1._dp - GzC%eps)*(log(C%o2)+ln106) ) )
      ! Now add in the bits that are reaction specific
      GzR%dG(1) = dG_ri + rV%RkJ*rV%T_K*(C%C_Phy/C%PhyS*(log(C%C_L)+ln106) - ( (log(C%PhyS)+ln106) + C%C_Phy/C%PhyS*(log(C%C_Phy)+ln106) ) ) ! PhyS predation
      GzR%dG(2) = dG_ri + rV%RkJ*rV%T_K*(C%C_GSB/C%GSBS*(log(C%C_L)+ln106) - ( (log(C%GSBS)+ln106) + C%C_GSB/C%GSBS*(log(C%C_GSB)+ln106) ) ) ! GSBS predation
      GzR%dG(3) = dG_ri - rV%RkJ*rV%T_K*(log(C%GzS) +ln106)  ! GzS  canibilization
      GzR%dG(4) = dG_ri - rV%RkJ*rV%T_K*(log(C%AGzS)+ln106)  ! AGzS predation
      GzR%dG(5) = dG_ri - rV%RkJ*rV%T_K*(log(C%BacS)+ln106)  ! BacS predation
      GzR%dG(6) = dG_ri - rV%RkJ*rV%T_K*(log(C%SRBS)+ln106)  ! SRBS predation
      GzR%dG(7) = dG_ri - rV%RkJ*rV%T_K*(log(C%PhS) +ln106)  ! PhS  predation
      GzR%dG(8) = dG_ri - rV%RkJ*rV%T_K*(log(C%SOxS)+ln106)  ! SOx  predation
      GzR%nei = 4._dp*GzR%aCi ! Electron transfer same for all reactions
      ! Reaction rates
      do i=1,8; GzR%F_T(i) = F_Thermo(GzR%dG(i), GzR%nei, rV%T_K); end do ! Thermo driver
      ! kinetic drivers
      F_o2 = C%o2/(C%o2 + rV%kappa*GzC%eps**4)
      GzR%F_K(1) = F_o2 * C%PhyS/(C%PhyS + rV%kappa*GzC%eps**4)
      GzR%F_K(2) = F_o2 * C%GSBS/(C%GSBS + rV%kappa*GzC%eps**4)
      GzR%F_K(3) = F_o2 * C%GzS /(C%GzS  + rV%kappa*GzC%eps**4)
      GzR%F_K(4) = F_o2 * C%AGzS/(C%AGzS + rV%kappa*GzC%eps**4)
      GzR%F_K(5) = F_o2 * C%BacS/(C%BacS + rV%kappa*GzC%eps**4)
      GzR%F_K(6) = F_o2 * C%SRBS/(C%SRBS + rV%kappa*GzC%eps**4)
      GzR%F_K(7) = F_o2 * C%PhS /(C%PhS  + rV%kappa*GzC%eps**4)
      GzR%F_K(8) = F_o2 * C%SOxS/(C%SOxS + rV%kappa*GzC%eps**4)

      ri = rV%nuStar*GzC%eps**2*C%GzS
      GzR%r(1) = ri*GzC%omg(1)*GzR%F_K(1) * GzR%F_T(1)
      GzR%r(2) = ri*GzC%omg(2)*GzR%F_K(2) * GzR%F_T(2)
      GzR%r(3) = ri*GzC%omg(3)*GzR%F_K(3) * GzR%F_T(3)
      GzR%r(4) = ri*GzC%omg(4)*GzR%F_K(4) * GzR%F_T(4)
      GzR%r(5) = ri*GzC%omg(5)*GzR%F_K(5) * GzR%F_T(5)
      GzR%r(6) = ri*GzC%omg(6)*GzR%F_K(6) * GzR%F_T(6)
      GzR%r(7) = ri*GzC%omg(7)*GzR%F_K(7) * GzR%F_T(7)
      GzR%r(8) = ri*GzC%omg(8)*GzR%F_K(8) * GzR%F_T(8)
      ! Entropy Production
      GzR%sig(1:8) = - GzR%r(1:8)*GzR%dG(1:8)*rV%csAtx/rV%T_K
      return
   end subroutine Gz_BioS
   
   subroutine AGz_bioS (t, x, C, AGzC, rV, AGzR)
      ! This routine is use to model anaerobic grazers (4: AGz)
      ! See Notes for Frontiers MS
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a, etc) and biological Structure conc (conc. in mmol/m3)
      type(AGzCntl)       , intent(in) :: AGzC  ! Control variables, eps, omg(1:7) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(AGzRxns)       , intent(out):: AGzR  ! all outputs associated with the anaerobic grazer reactions.  See module typeVars
      
      ! local declarations
      integer i
      real(dp) dG0_ri, dG_ri ! Free energy for reaction i, which is the consumption rate of BioS(i) (kJ/mol or J/mmol)
      real(dp) ri ! for most of the reaction rates
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
      real(dp) F_so4
      
      ! Catabolic stiochiometric coef.
      AGzR%aCi = (  4._dp + rV%alp - 2._dp*rV%bet - 3._dp*rV%gam + 5._dp*rV%del - 4._dp*AGzC%eps)/8._dp
      AGzR%bCi = (- 2._dp + rV%alp                - 3._dp*rV%gam - 3._dp*rV%del                 )/2._dp
      
      ! Free energy of reactions (Note, this does not depend on C_Phy or C_GSB), so these are all the same
      dG0_ri = AGzC%eps*rV%dGf0_bioS + (1._dp - AGzC%eps)**2*rV%dGf0_dic + AGzC%eps*(1._dp - AGzC%eps)*rV%dGf0_ch2o + rv%del*(1._dp - AGzC%eps)*rV%dGf0_h3po4 &
             & + rv%gam*(1._dp - AGzC%eps)*rV%dGf0_nh3 + AGzR%bCi*(1._dp - AGzC%eps)*rV%dGf0_h2o + AGzR%aCi*(1._dp - AGzC%eps)*rV%dGf0_h2saq &
             & - ( rV%dGf0_bioS + AGzR%aCi*(1._dp - AGzC%eps)*rV%dGf0_h2so4 )
      ! First calculate the quantity that is the same for all reactions.
      dG_ri  = dG0_ri + rV%RkJ*rV%T_K*( AGzC%eps*(log(C%AGzS)+ln106) + (1._dp - AGzC%eps)**2*(log(C%dic)+ln106) + AGzC%eps*(1._dp - AGzC%eps)*(log(C%C_D)+ln106) &
             + rv%del*(1._dp - AGzC%eps)*((1._dp - AGzC%eps)*(log(C%po4)+ln106) + AGzC%eps*(log(C%P_D)+ln106)) &
             + rv%gam*(1._dp - AGzC%eps)*((1._dp - AGzC%eps)*(log(C%nh3)+ln106) + AGzC%eps*(log(C%N_D)+ln106)) &
             + AGzR%aCi*(1._dp - AGzC%eps)*(log(C%h2s)+ln106) - (AGzR%aCi*(1._dp - AGzC%eps)*(log(C%so4)+ln106) ) )
      ! Now add in the bits that are reaction specific
      AGzR%dG(1) = dG_ri + rV%RkJ*rV%T_K*( C%C_Phy/C%PhyS*(log(C%C_L)+ln106) - ( (log(C%PhyS)+ln106) + C%C_Phy/C%PhyS*(log(C%C_Phy)+ln106) ) ) ! PhyS predation
      AGzR%dG(2) = dG_ri + rV%RkJ*rV%T_K*( C%C_GSB/C%GSBS*(log(C%C_L)+ln106) - ( (log(C%GSBS)+ln106) + C%C_GSB/C%GSBS*(log(C%C_GSB)+ln106) ) )! GSBS predation
      AGzR%dG(3) = dG_ri - rV%RkJ*rV%T_K*(log(C%GzS) +ln106)  ! GzS  predation
      AGzR%dG(4) = dG_ri - rV%RkJ*rV%T_K*(log(C%AGzS)+ln106)  ! AGzS canibilization
      AGzR%dG(5) = dG_ri - rV%RkJ*rV%T_K*(log(C%BacS)+ln106)  ! BacS predation
      AGzR%dG(6) = dG_ri - rV%RkJ*rV%T_K*(log(C%SRBS)+ln106)  ! SRBS predation
      AGzR%dG(7) = dG_ri - rV%RkJ*rV%T_K*(log(C%PhS) +ln106)  ! PhS  predation
      AGzR%dG(8) = dG_ri - rV%RkJ*rV%T_K*(log(C%SOxS)+ln106)  ! SOx  predation
      AGzR%nei = 8._dp*AGzR%aCi ! Electron transfer same for all reactions
      ! Reaction rates
      do i=1,8; AGzR%F_T(i) = F_Thermo(AGzR%dG(i), AGzR%nei, rV%T_K); end do ! Thermo driver
      ! kinetic drivers
      F_so4 = C%so4/(C%so4 + rV%kappa*AGzC%eps**4)
      AGzR%F_K(1) = F_so4 * C%PhyS/(C%PhyS + rV%kappa*AGzC%eps**4)
      AGzR%F_K(2) = F_so4 * C%GSBS/(C%GSBS + rV%kappa*AGzC%eps**4)
      AGzR%F_K(3) = F_so4 * C%GzS /(C%GzS  + rV%kappa*AGzC%eps**4)
      AGzR%F_K(4) = F_so4 * C%AGzS/(C%AGzS + rV%kappa*AGzC%eps**4)
      AGzR%F_K(5) = F_so4 * C%BacS/(C%BacS + rV%kappa*AGzC%eps**4)
      AGzR%F_K(6) = F_so4 * C%SRBS/(C%SRBS + rV%kappa*AGzC%eps**4)
      AGzR%F_K(7) = F_so4 * C%PhS /(C%PhS  + rV%kappa*AGzC%eps**4)
      AGzR%F_K(8) = F_so4 * C%SOxS/(C%SOxS + rV%kappa*AGzC%eps**4)
         
      ri = rV%nuStar*AGzC%eps**2*C%AGzS
      AGzR%r(1) = ri*AGzC%omg(1)*AGzR%F_K(1) * AGzR%F_T(1)
      AGzR%r(2) = ri*AGzC%omg(2)*AGzR%F_K(2) * AGzR%F_T(2)
      AGzR%r(3) = ri*AGzC%omg(3)*AGzR%F_K(3) * AGzR%F_T(3)
      AGzR%r(4) = ri*AGzC%omg(4)*AGzR%F_K(4) * AGzR%F_T(4)
      AGzR%r(5) = ri*AGzC%omg(5)*AGzR%F_K(5) * AGzR%F_T(5)
      AGzR%r(6) = ri*AGzC%omg(6)*AGzR%F_K(6) * AGzR%F_T(6)
      AGzR%r(7) = ri*AGzC%omg(7)*AGzR%F_K(7) * AGzR%F_T(7)
      AGzR%r(8) = ri*AGzC%omg(8)*AGzR%F_K(8) * AGzR%F_T(8)      
      
      ! Entropy Production
      AGzR%sig(1:8) = - AGzR%r(1:8)*AGzR%dG(1:8)*rV%csAtx/rV%T_K
      return
   end subroutine AGz_BioS

   subroutine Bac_bioS (t, x, C, BacC, rV, BacR)
      ! This routine is use to model heterotrophic bacteria (5: Bac)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! Three reactions are considered here: 1) Bacterial growth, 2) C_D decomposition, 3) P_D decomposition
      use realPrec
      use typeVars
      use functions, only: F_Thermo, stepUp
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a) and biological Structure conc (conc. in mmol/m3)
      type(BacCntl)       , intent(in) :: BacC  ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(BacRxns)       , intent(out):: BacR  ! all outputs associated with the Bac reactions.  See module typeVars
      
      ! local declarations
      real(dp) dG0_rA1, dG_rA1, dG0_rC1, dG_rC1 ! free energy for anabolic and catabolic reactions for r1
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Calculate heterogrophic growth, r1
      ! C_L + eps(gam NH3 + del H2PO4) + (1-eps)o2 -> eps aA1 Bac + eps bA1 H2O + [2 - eps(aA1-1)]dic
      ! calculate reaction coefficients, but only anabolic reactions have them
      
      BacR%aA1 = (4._dp + 3._dp*rV%gam - 5._dp*rV%del)/(4._dp + rV%alp - 2._dp*rV%bet)
      BacR%bA1 = (4._dp - 2._dp*rV%alp + 9._dp*rV%gam - 3._dp*rV%bet*rV%gam + rV%del + 4._dp*rV%alp*rV%del - 3._dp*rV%bet*rV%del)/(4._dp + rV%alp - 2._dp*rV%bet)
      ! anabolic reaction thermodynamics
      dG0_rA1 = (BacR%aA1*rV%dGf0_bioS + (1._dp-BacR%aA1)*rV%dGf0_dic + BacR%bA1*rV%dGf0_h2o) - (rV%dGf0_ch2o + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4)
      dG_rA1  = dG0_rA1 + rV%RkJ*rV%T_K*( BacR%aA1*(log(C%BacS)+ln106) + (1._dp-BacR%aA1)*(log(C%dic)+ln106) - (log(C%C_L)+ln106) - rV%gam*(log(C%nh3)+ln106) - rV%del*(log(C%po4)+ln106) )

      ! catabolic reaction thermodynamics
      dG0_rC1 = rV%dGf0_dic - ( rV%dGf0_ch2o + rV%dGf0_o2aq )
      dG_rC1  = dG0_rC1 + rV%RkJ*rV%T_K*( (log(C%dic)+ln106) - (log(C%C_L)+ln106) - (log(C%o2)+ln106) )
      ! Free energy of combined whole reaction
      BacR%dG(1) = BacC%eps*dG_rA1 + (1._dp - BacC%eps)*dG_rC1
      ! Reaction (1) rate
      BacR%ne = 4. ! electrons transfered in catabolic reaction
      BacR%F_T(1) = F_Thermo(BacR%dG(1), BacR%ne, rV%T_K)
      BacR%F_K(1) = ( C%C_L/(C%C_L + rV%kappa*BacC%eps**4) )*( C%po4/rV%del/(C%po4/rV%del + rV%kappa*BacC%eps**4)  )*( C%o2/(C%o2 + rV%kappa*BacC%eps**4)  ) ! don't include NH3 limatations yet
      BacR%r(1) = rV%nuStar*BacC%eps**2*BacC%omg(1)*C%BacS*BacR%F_T(1)*BacR%F_K(1)
      ! Entropy production
      BacR%sig(1) = - BacR%r(1)*BacR%dG(1)*rV%csAtx/rV%T_K
      
      ! Reaction 2, C_D decompositoin to C_L.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! C_D -> C_L
      BacR%dG(2) = rV%RkJ*rV%T_K*log(C%C_L/C%C_D) ! This will be negligable, but easy to calculate so keep
      BacR%F_K(2) = C%C_D/(C%C_D + rV%kappa*BacC%eps**4)
      BacR%F_T(2) = stepUp(-BacR%dG(2), 0.3_dp, 20._dp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      BacR%r(2) = rV%nuDet*BacC%eps**2*BacC%omg(2)*C%BacS*BacR%F_T(2)*BacR%F_K(2)
      BacR%sig(2) = - BacR%r(2)*BacR%dG(2)*rV%csAtx/rV%T_K
      
      ! Reaction 3, P_D decompositoin to H2PO4.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! P_D -> H3PO4
      BacR%dG(3) = rV%RkJ*rV%T_K*log(C%po4/C%P_D) ! This will be negligable, but easy to calculate so keep
      BacR%F_K(3) = C%P_D/(C%P_D + rV%kappa*BacC%eps**4)
      BacR%F_T(3) = stepUp(-BacR%dG(3), 0.3_dp, 20._dp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      BacR%r(3) = rV%nuDet*BacC%eps**2*BacC%omg(3)*C%BacS*BacR%F_T(3)*BacR%F_K(3)
      BacR%sig(3) = - BacR%r(3)*BacR%dG(3)*rV%csAtx/rV%T_K
      
      return
   end subroutine Bac_bioS
   
   subroutine SRB_bioS (t, x, C, SRBC, rV, SRBR)
      ! This routine is use to model sulfate reducing bacteria (6: SRB)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! Three reactions are considered here: 1) Bacterial growth, 2) C_D decomposition, 3) P_D decomposition
      use realPrec
      use typeVars
      use functions, only: F_Thermo, stepUp
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a) and biological Structure conc (conc. in mmol/m3)
      type(SRBCntl)       , intent(in) :: SRBC  ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(SRBRxns)       , intent(out):: SRBR  ! all outputs associated with the SRB reactions.  See module typeVars
      
      ! local declarations
      real(dp) dG0_rA1, dG_rA1, dG0_rC1, dG_rC1 ! free energy for anabolic and catabolic reactions for r1
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Calculate SRB growth, r1
      ! (1+eps n)C_L + eps(gam NH3 + del H2PO4) + [1/2 + eps(aA +1/2n-1/2)]h2so4 -> eps SRB + [1+eps(n-1)]dic + [1/2 + eps(aA +1/2n-1/2)]h2s + eps bA h2o
      ! calculate reaction coefficients, but only anabolic reactions have them
      SRBR%aA1 = (      - rV%alp + 2._dp*rV%bet + 3._dp*rV%gam - 5._dp*rV%del)/8._dp
      SRBR%bA1 = (2._dp - rV%alp                + 3._dp*rV%gam + 3._dp*rV%del)/2._dp
      ! anabolic reaction thermodynamics
      dG0_rA1 = (rV%dGf0_bioS +  SRBR%bA1*rV%dGf0_h2o + SRBR%aA1*rV%dGf0_h2saq) - (rV%dGf0_ch2o + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4 + SRBR%aA1*rV%dGf0_h2so4)
      dG_rA1  = dG0_rA1 + rV%RkJ*rV%T_K*( (log(C%SRBS)+ln106) + SRBR%aA1*(log(C%h2s)+ln106) - (log(C%C_L)+ln106) - rV%gam*(log(C%nh3)+ln106) &
              & - rV%del*(log(C%po4)+ln106) - SRBR%aA1*(log(C%so4)+ln106)  )
      ! catabolic reaction thermodynamics
      dG0_rC1 = rV%dGf0_dic + rV%dGf0_h2saq/2._dp - ( rV%dGf0_ch2o + rV%dGf0_h2so4/2._dp )
      dG_rC1 = dG0_rC1 + rV%RkJ*rV%T_K*( (log(C%dic)+ln106) + (log(C%h2s)+ln106)/2._dp - (log(C%C_L)+ln106) - (log(C%so4)+ln106)/2._dp )
      ! coupling of anabolic to catabolic reaction
      SRBR%n1 = - dG_rA1/dG_rC1
      ! Free energy of combined whole reaction
      SRBR%dG(1) = (1._dp - SRBC%eps)*dG_rC1
      ! Reaction (1) rate
      SRBR%ne = 4. ! electrons transfered in catabolic reaction
      SRBR%F_T(1) = F_Thermo(SRBR%dG(1), SRBR%ne, rV%T_K)
      SRBR%F_K(1) = ( C%C_L/(C%C_L + rV%kappa*SRBC%eps**4) )*( C%po4/rV%del/(C%po4/rV%del + rV%kappa*SRBC%eps**4)  )*( C%so4/(C%so4 + rV%kappa*SRBC%eps**4)  ) ! don't include NH3 limatations yet
      SRBR%r(1) = rV%nuStar*SRBC%eps**2*SRBC%omg(1)*C%SRBS*SRBR%F_T(1)*SRBR%F_K(1)
      ! Entropy production
      SRBR%sig(1) = - SRBR%r(1)*SRBR%dG(1)*rV%csAtx/rV%T_K
      
      ! Reaction 2, C_D decompositoin to C_L.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! C_D -> C_L
      SRBR%dG(2) = rV%RkJ*rV%T_K*log(C%C_L/C%C_D) ! This will be negligable, but easy to calculate so keep
      SRBR%F_K(2) = C%C_D/(C%C_D + rV%kappa*SRBC%eps**4)
      SRBR%F_T(2) = stepUp(-SRBR%dG(2), 0.3_dp, 20._dp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      SRBR%r(2) = rV%nuDet*SRBC%eps**2*SRBC%omg(2)*C%SRBS*SRBR%F_T(2)*SRBR%F_K(2)
      SRBR%sig(2) = - SRBR%r(2)*SRBR%dG(2)*rV%csAtx/rV%T_K
      
      ! Reaction 3, P_D decompositoin to H2PO4.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! P_D -> H3PO4
      SRBR%dG(3) = rV%RkJ*rV%T_K*log(C%po4/C%P_D) ! This will be negligable, but easy to calculate so keep
      SRBR%F_K(3) = C%P_D/(C%P_D + rV%kappa*SRBC%eps**4)
      SRBR%F_T(3) = stepUp(-SRBR%dG(3), 0.3_dp, 20._dp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      SRBR%r(3) = rV%nuDet*SRBC%eps**2*SRBC%omg(3)*C%SRBS*SRBR%F_T(3)*SRBR%F_K(3)
      SRBR%sig(3) = - SRBR%r(3)*SRBR%dG(3)*rV%csAtx/rV%T_K
      
      return
   end subroutine SRB_bioS
   
   subroutine Ph_bioS (t, x, C, PhC, rV, PhR)
      ! This routine is use to model Photoheterotrophs (7: Ph)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! Only reactions is considered here: Conversion of C_L to PhS
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables (o2 , dic , po4 , C_a) and biological Structure conc (conc. in mmol/m3)
      type(PhCntl )       , intent(in) :: PhC   ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(PhRxns )       , intent(out):: PhR   ! all outputs associated with the Ph reactions.  See module typeVars
      
      ! local declarations
      real(dp) dG0_rA1, dG_rA1 ! free energy for anabolic and catabolic reactions for r1
      real(dp) delI_Ph
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Calculate photoheterogrophic growth, r1
      ! eps(C_L + eps(gam NH3 + del H2PO4) + hv -> eps( aA1 Ph + bA1 H2O + (1 - aA1-1) )dic
      ! calculate reaction coefficients, but only anabolic reactions have them 
      PhR%aA1 = (4._dp + 3._dp*rV%gam - 5._dp*rV%del)/(4._dp + rV%alp - 2._dp*rV%bet)
      PhR%bA1 = (4._dp - 2._dp*rV%alp + 9._dp*rV%gam - 3._dp*rV%bet*rV%gam + rV%del + 4._dp*rV%alp*rV%del - 3._dp*rV%bet*rV%del)/(4._dp + rV%alp - 2._dp*rV%bet)
      ! anabolic reaction thermodynamics
      dG0_rA1 = (PhR%aA1*rV%dGf0_bioS + (1._dp-PhR%aA1)*rV%dGf0_dic + PhR%bA1*rV%dGf0_h2o) - (rV%dGf0_ch2o + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4)
      dG_rA1  = dG0_rA1 + rV%RkJ*rV%T_K*( PhR%aA1*(log(C%PhS)+ln106) + (1._dp-PhR%aA1)*(log(C%dic)+ln106) - (log(C%C_L)+ln106) - rV%gam*(log(C%nh3)+ln106) - rV%del*(log(C%po4)+ln106) )
      ! Light capture
      delI_Ph  = rV%k_p*C%PhS*rV%Itx ! note, assume all PhS allocated to light harvesting for chemosynthesis.
      ! Free energy of combined whole reaction
      PhR%dG(1) = PhC%eps*dG_rA1 + rV%dGr_Ggamma
      ! Reaction (1) rate
      PhR%ne = 1. ! electrons transfered by photons (see pg 47 of notes)
      PhR%F_T(1) = F_Thermo(PhR%dG(1), PhR%ne, rV%T_K)
      PhR%F_K(1) = ( C%C_L/(C%C_L + rV%kappa*PhC%eps**4) )*( C%po4/rV%del/(C%po4/rV%del + rV%kappa*PhC%eps**4) ) ! don't include NH3 limatations yet
      PhR%r(1) = delI_Ph*PhR%F_T(1)*PhR%F_K(1)
      ! Entropy production
      PhR%sig(1) = - ( delI_Ph*rV%dGr_Ggamma + PhC%eps*PhR%r(1)*dG_rA1 )*rV%csAtx/rV%T_K
      
      return
   end subroutine Ph_bioS

   subroutine SOx_bioS (t, x, C, SOxC, rV, SOxR)
      ! This routine is use to model sulfur bacteria (8: SOx)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS, pp. 85-87
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      use realPrec
      use typeVars
      use functions, only: F_Thermo
      implicit none
      real(dp)            , intent(in) :: t     ! time (d)
      real(dp)            , intent(in) :: x     ! location (m)
      type(concentrations), intent(in) :: C     ! conc variables and biological Structure conc (conc. in mmol/m3)
      type(SOxCntl )      , intent(in) :: SOxC  ! Control variables, eps, omg(1:2) (the capital omegas that sum to 1.0, not lower case omega, w)
      type(rxnVariables)  , intent(in) :: rV    ! Parameters needed by for calculations.  See rxnVariables type in module typeVars.
      type(SOxRxns )      , intent(out):: SOxR  ! all outputs associated with the Ph reactions.  See module typeVars
      
      ! local declarations
      real(dp) dG0_rA1, dG_rA1, dG0_rC1, dG_rC1 ! free energy for anabolic and catabolic reactions for r1
      real(dp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M
            
      ! Calculate SOx growth, r1
      ! eps H2CO3 + [eps(a+n-1)+1]H2S + 2[eps(n-1)+1]O2 + eps gam NH3 + eps del H3PO4 -> eps SOxS + [eps(a+n-1)+1]H2SO4 + b H2O
      ! calculate reaction coefficients, but only anabolic reactions have them 
      SOxR%aA1 = 0.1_dp*(4._dp +       rV%alp - 2._dp*rV%bet - 3._dp*rV%gam +  5._dp*rV%del)
      SOxR%bA1 = 0.2_dp*(7._dp - 2._dp*rV%alp -       rV%bet + 6._dp*rV%gam + 10._dp*rV%del)
      ! anabolic reaction thermodynamics
      dG0_rA1 = (rV%dGf0_bioS + SOxR%bA1*rV%dGf0_h2o + SOxR%aA1*rV%dGf0_h2so4) - (rV%dGf0_dic + SOxR%aA1*rV%dGf0_h2saq + rV%gam*rV%dGf0_nh3 + rV%del*rV%dGf0_h3po4)
      dG_rA1  = dG0_rA1 + rV%RkJ*rV%T_K*( (log(C%SOxS)+ln106) + SOxR%aA1*(log(C%so4)+ln106) - (log(C%dic)+ln106) - SOxR%aA1*(log(C%h2s)+ln106) - rV%gam*(log(C%nh3)+ln106) &
              & - rV%del*(log(C%po4)+ln106) )
      ! catabolic reaction
      dG0_rC1 = rV%dGf0_h2so4 - ( rV%dGf0_h2saq + 2._dp*rV%dGf0_o2aq )
      dG_rC1 = dG0_rC1 + rV%RkJ*rV%T_K*( (log(C%so4)+ln106) - (log(C%h2s)+ln106) - 2._dp*(log(C%o2)+ln106) )
      ! coupling of anabolic to catabolic reaction
      SOxR%n1 = - dG_rA1/dG_rC1
      ! Free energy of combined whole reaction      
      SOxR%dG(1) = (1._dp - SOxC%eps)*dG_rC1
      ! Reaction (1) rate
      SOxR%ne = 8. ! electrons transfered by photons (see pg 86 of notes)
      SOxR%F_T(1) = F_Thermo(SOxR%dG(1), SOxR%ne, rV%T_K)
      SOxR%F_K(1) = (C%dic/(C%dic + rV%kappa*SOxC%eps**4))*(C%h2s/(C%h2s + rV%kappa*SOxC%eps**4))*(C%o2/(C%o2 + rV%kappa*SOxC%eps**4))*(C%po4/rV%del/(C%po4/rV%del + rV%kappa*SOxC%eps**4)) ! don't include NH3 limatations yet
      SOxR%r(1) = rV%nuStar*SOxC%eps**2*C%SOxS*SOxR%F_T(1)*SOxR%F_K(1)
      ! Entropy production
      SOxR%sig(1) = - SOxR%r(1)*SOxR%dG(1)*rV%csAtx/rV%T_K    
      return
   end subroutine SOx_bioS
   
   Subroutine parSur (time,dlat,solarC,par)
      use realPrec
      Implicit None
      real(dp) time, dlat, par, solarC
      !
      !     This routine the light level at the surface of the earth given
      !     time and latitude (and some other stuff)
      !
      !     USAGE - Call parSur (time,dlat,solarC,par)
      !
      !     INPUT
      !        TIME     REAL*8. Time of year in days (Julian days)
      !        dlat     Real*8  Latitude in degrees
      !        solarC   Real*8  Solar constant, in whatever units you want par in.
      !                         e.g.  1353. W/m2
      !                               6216. microEinstein/s/m2 Note, 0.75 of this seems to give more reasonable values (2000 vs 2650 on day 177.5)
      !
      !     OUTPUT
      !        par      REAL*8. photosynthetic available radiation.
      !
      !     MODIFIED:   31Jan97, 7Mar2016

      !     Local Declarations
      real(dp) pi, rlat, hour, dec, ha, cosz, aa           
      !real(dp) sigma ! this is an exponential weighting coef. to allow a smooth transition to night.
      !      
      !     Declination and cos of Zenith angle.
      !     Ref: Brock (1981) Ecol.Model. 14:1-19, changed to radians
      pi = 3.141592654D0
      rlat = dlat*pi/180.D0
      hour = (time - aint(time))*24.D0
      dec = 0.4093D0*sin(2.*pi*(284.D0+time)/365.D0)
      ha = (hour - 12.D0)*pi/12.D0                                           
      cosz = Dsin(dec)*Dsin(rlat) + Dcos(dec)*Dcos(rlat)*Dcos(ha)           
      !     See if it is night time    
      !If (cosz <= 0._dp) cosz = 0.0D0

      ! see if exponential weighting that removes the discontinuity helps numerical solution
      ! see C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\Mma\SmoothPARfcn.nb
      !sigma = 50._dp
      !If (cosz <= 0._dp) then
      !   cosz = 0.0D0
      !else
      !   cosz = cosz*( 1._dp - exp(-sigma*cosz) )/( 1._dp + exp(-sigma*cosz) )
      !end if            
      cosz = max(0._dp, cosz) ! the above exponential function did not seem to matter, so just use this
      
      ! PAR to Total radiation at surface: a = Epar(0+)/Etot(0+)
      ! Ref: Baker & Frouin (1987) L&O 32:1370-1377     
      aa = 0.45 
      par = solarC*aa*cosz      
      return
   end  
   
   subroutine interp1D (t, ndim_t, tvec, nlow, delt)
      ! This routine determines nlow such that:
      !    tvec(nlow) <= t <= tvec(nlow+1)
      ! Note, if t = tvec(ndim_t), then nlow is set to ndim_t-1 
      ! delt is given by:
      !    delt = [t - tvec(nlow)]/[tvec(nlow+1)-tvec(nlow)]
      ! that can be used for linear interpolation

      ! Input
      !    t        time that is to be found in tvec
      !    ndim_t   Number of used elements of tvec
      !    tvec     Vector containing increasing values of t
      ! Output
      !    nlow     Value such that: tvec(nlow) <= t <= tvec(nlow+1)
      !    delt     value such that: [t - tvec(nlow)]/[tvec(nlow+1)-tvec(nlow)]
      use realPrec
      Implicit None              
      Integer ndim_t, nlow
      real(dp) t, tvec(ndim_t), delt

      ! Local Declarations:
      Integer max_itr, nhi, nmid, i

      ! If t is outside tvec range, just return upper or lower bound
      If (t <= tvec(1)) Then
         nlow = 1
         delt = 0.0
         Return
      End If
      If (t >= tvec(ndim_t)) Then
         nlow = ndim_t-1
         delt = 1.0
         Return
      End If

      ! find nlow by binary search.
      max_itr = Int( Log(real(ndim_t))/Log(2.0_dp) + 2.0_dp )
      nlow = 1
      nhi = ndim_t      
      nmid = (nhi + nlow)/2 
      Do i=1,max_itr
         If (nhi == nlow + 1) then
            delt = (t - tvec(nlow))/( tvec(nlow+1)-tvec(nlow) )
            return
         end if         
         If (t < tvec(nmid)) Then
            nhi = nmid
         Else
            nlow = nmid
         End If 
         nmid = (nhi + nlow)/2         
      End Do                     
      Write(*,*) 'interp1D::Search failed (this should not happen)'  
      Return
   end subroutine interp1D

   subroutine interp2D (nt, nx, nk, tgrid, xgrid, utxk, t, x, u)
      ! This routine interpolates in 2D (t,x) in utxk.  In this version, utxk contains
      ! a vector of variables, k, where each variable spaces a 2D space (t,x).
      ! This space is a regular grid, given by the vectors tgrid and xgrid
      ! Given a value of t and x, this routine returns a vector, u(nk), that is interploated to (t,x) 
      use realPrec
      implicit none
      integer,  intent(in) :: nt, nx, nk     ! demenstion of utxk array
      real(dp), intent(in) :: tgrid(nt)      ! grid points in t
      real(dp), intent(in) :: xgrid(nx)      ! grid points in x
      real(dp), intent(in) :: utxk(nt,nx,nk) ! array of u values
      real(dp), intent(in) :: t              ! time to interploate to
      real(dp), intent(in) :: x              ! space to interploate to
      real(dp), intent(out):: u(nk)          ! values of u for all k at (t,x)
    
      ! Local declaratiosn
      integer ntlow, nxlow, k
      real(dp) delt, delx, uxti, uxtip1
    
      ! First get interpolated locations in t and x dimensions
      call interp1D (t, nt, tgrid, ntlow, delt)
      call interp1D (x, nx, xgrid, nxlow, delx)
    
      do k=1, nk
         uxti   = utxk(ntlow  ,nxlow, k) + delx*(utxk(ntlow  ,nxlow+1, k) - utxk(ntlow  ,nxlow, k))
         uxtip1 = utxk(ntlow+1,nxlow, k) + delx*(utxk(ntlow+1,nxlow+1, k) - utxk(ntlow+1,nxlow, k))
         u(k) = uxti + delt*(uxtip1 - uxti)
      end do
      return
   end subroutine interp2D
   
   subroutine spaceNodes ()
      ! this routine spaces the splne nodes exponentially; however, if the
      ! sigma weighting function is close to linear, nodes are spaced evently.
      ! NOTE tnode(0) must be set prior to calling this routine, and it should equal t at start of interval 
      use parameters
      implicit none
   
      integer i
      real(8) wS
   
      if (wInfinity >= 0.999) then
         ! just space the knots equally.
         do i=1, ncntl(1)
            tOCnodes(i) = tOCnodes(i-1) + deltInfinity/real(ncntl(1))
         end do      
         return
      end if
      ! Space knots exponentially (see function weightS, as this is its inverse)
      wS = 1.0d0
      do i=1,ncntl(1)-1
         wS = wS - (1.d0 - wInfinity)/real(ncntl(1))
         tOCnodes(i) = deltInfinity*log(wS)/log(wInfinity) + tOCnodes(0)
      end do 
      tOCnodes(ncntl(1)) = tOCnodes(0) + deltInfinity !just avoids any rounding errors in above expression  
      return
   end subroutine spaceNodes    
   
   subroutine omega_free (n, w, omg)
      ! this routine calclates omega (capital) from the lower case omega based on sorting
      ! of w then calculating distances between w_i.  
      ! Note, if there is only 1 sub reaction, then omg = w = 1, so there is no need
      ! to call this routine, but it is robust to such a call.
      use realPrec
      use interfaces
      use typeVars, only: maxSubRxn
      implicit none
      integer, intent(in)  :: n              ! number of sub reactions for a biological structure
      !real(dp), intent(in) :: w(max(1,n-1))  ! lower case omega, control variables
      !real(dp), intent(out):: omg(n)         ! upper case omega, these sum to 1, which is where the nth omega comes from
      real(dp), intent(in) :: w(*)  ! lower case omega, control variables
      real(dp), intent(out):: omg(*)         ! upper case omega, these sum to 1, which is where the nth omega comes from
      ! Local declaration
      integer i, l
      real(dp) w1(maxSubRxn+1) ! By using a defined size vector prevents it from being automatic and constantly allocated/deallocated 

      omg(1) = 1.0_dp ! set this in case n == 1
      if (n == 1) return
      ! For some reason this construction took a lot of CPU time, or at least that's what it looks like
      !w1(1:n+1) = [0._dp, w(1:n-1), 1._dp] ! add 0 and 1 to ends
      
      ! note, if the lower bound on omega is greater than zero, it is possible for some omega to be zero.  However, this may
      ! not necessarily be a problem, so check into it.
      ! omg-1 will always be >= lowerBnd if w >= lowerBnd
      ! w(1) should be set to lowerBnd
      ! w(n+1) = upperBnd  
      ! should omg be scaled so that they do very between 0 and 1?
      ! Currently, upperBnd is 1, and lowerBnd is 10^-8, so it's not really a problem.

      w1(1) = 0_dp
      w1(2:n) = w(1:n-1)
      w1(n+1) = 1_dp
      call insertion_sort(n, w1) ! sort w1
      omg(1:n) = w1(2:n+1) - w1(1:n)
      return
   end subroutine omega_free
   
   pure subroutine insertion_sort(n, a)
      ! uses insertion sorting
      ! from here https://rosettacode.org/wiki/Category:Sorting_Algorithms      
      ! but also see https://en.wikipedia.org/wiki/Sorting_algorithm
      use realPrec
      implicit none
      integer, intent(in)                    :: n
      real(dp), intent(in out), dimension(:) :: a
      real(dp) :: temp
      integer :: i, j
 
      do i = 2, n
         j = i - 1
         temp = a(i)
         do while (j>=1 .and. a(j)>temp) ! Strange, if I remove j>=1, it runs slower
            a(j+1) = a(j)
            j = j - 1
            if (j<=0) exit ! this prevents do while test of a(0)
         end do
         a(j+1) = temp
      end do
   end subroutine insertion_sort     
   
   !-----------------------------------------------------------------------
   ! Routines below are provide the definition of the PDE for BACOLI95
   !-----------------------------------------------------------------------
   subroutine PDEs(t, xd, u, ux, uxx, fval, npde)
      ! This routine defines the system of PDEs
      ! dc/dt = F(t, x, c, cx, cxx), were the function is passed back in fval
      !
      use realPrec
      use parameters
      use functions
      use typeVars
      implicit none
      integer, intent(in)     ::  npde        ! number of PDEs
      real(dp), intent(in)    ::  t           ! time (d)
      real(dp), intent(inout) ::  xd          ! spatial coord non-dimensional
      real(dp), intent(inout) ::  u(npde)     ! value of c at (t,x) (non-dimensional)
      real(dp), intent(inout) ::  ux(npde)    ! first spatial derivative (dc/dx) (non-dimensional)
      real(dp), intent(inout) ::  uxx(npde)   ! second spatial derivative (d2c/dx2) (non-dimensional)
      real(dp), intent(out)   ::  fval(npde)  ! value of dc/dt = fval(t, x, c, cx, cxx) (units of 1/d)
      ! local declarations
      real(dp) sumSigmaDot(3) ! reaction entropy production for the three components (rxns, water, particles)
      real(dp) cntlTX(ncntl(3))
      real(dp) transport(maxPDE), D, dDdx, a, dadx, q, dqdx, qL, cL(maxPDE)
      type(concentrations) Cn
      type(rxnNet) rN    ! This setups all outputs from rxnProperties.  See typeVars module
      real(dp) x, c(maxPDE), cx(maxPDE), cxx(maxPDE) ! These automatic arrays may be put on the stack.  Use allocatable to put them on the heap.
      real(dp) vS(maxPDE) ! used to improve PDE startup/initialization. Spins up vSink
      real(dp) spinUp
      
      ! Convert dimensionales quantities to their dimensional forms
      x = xd*Xc
      c = u*Cc(1:npde)
      cx = ux*Cc(1:npde)/Xc
      cxx = uxx*Cc(1:npde)/(Xc*Xc)
      where (c < 0._dp) c = absZero ! prevent negative terms
      
      ! set names, and insure none are <= to absZero, and also makes it easy to refer to them
      call setCnames (c(1:nconc), c(nconc+1:nconc+nbioS), Cn)            
   
      ! get transport terms.  Note, all of these return fully dimensional quantities (NO normalization)
      call dispersion (t, x, D, dDdx) ! Dispersion and it's derivative (m2/d)
      call ccArea     (t, x, a, dadx) ! crossectional area (m2)
      call volFlow    (t, x, q, dqdx, qL, cL) ! volumetric flow rate (m3/d) and laterial inputs if any. qL is in m2/d
      spinUp = stepUp(t, tStep, tSig) ! This allows easier starting
      q = q*spinUp
      dqdx = spinUp*dqdx
      qL = qL*spinUp     
      vS = vSink*spinUp ! this slowly spins up sinking.

      
      ! calclate transport of constituent assocaite with advection and dispersion, full dimensions
      transport(1:npde) = D*cxx(1:npde) + (d*dadx/a + dDdx - q/a )*cx(1:npde) &
                        - vS(1:npde)*cx(1:npde) - dqdx/a*c(1:npde) - vS(1:npde)*c(1:npde)/a*dadx + qL/a*cL(1:npde)   

      ! Interpolate time in tOCnodes and x in xOCnodes
      call  interp2D (ncntl(1)+1, ncntl(2), ncntl(3), tOCnodes, xOCnodes, cntl, t, x, cntlTX)

      ! Calculate free energy of reactions, reaction rates, and thermodynamics for each catalyst
      call rxnProperties(t, x, c(1:nconc), c(nconc+1:nconc+nbioS), ncntl(3), cntlTX, rN, sumSigmaDot)
      
      ! Define PDE equations 
      ! salinity, sal
      fval(1) = transport(1) ! no reaction term on salinity
      ! dissolved oxygen, o2
      fval(2) = transport(2) + rN%PhyC%eps*rN%PhyR%r(1) - ( 1._dp + rN%PhyC%eps*(rN%PhyR%aA2 + rN%PhyR%n(2)-1._dp) )*rN%PhyR%r(2)  &
              - rN%GzR%aCi*(1._dp - rN%GzC%eps)*sum(rN%GzR%r(1:8)) - (1._dp - rN%BacC%eps)*rN%BacR%r(1) - 2._dp*(rN%SOxC%eps*(rN%SOxR%n1 - 1._dp) + 1._dp)*rN%SOxR%r(1) 
      ! dissolved inorganic carbon, dic
      fval(3) = transport(3) - rN%PhyC%eps*rN%PhyR%r(1) + (1._dp + rN%PhyC%eps*(rN%PhyR%n(2) - 1._dp))*rN%PhyR%r(2) - rN%GSBC%eps*rN%GSBR%r(1) &
              + (1._dp + rN%GSBC%eps*(rN%GSBR%n(2) - 1._dp))*rN%GSBR%r(2) + (1._dp - rN%GzC%eps)**2*sum(rN%GzR%r(1:8)) &
              + (1._dp - rN%AGzC%eps)**2*sum(rN%AGzR%r(1:8)) + (2._dp - rN%BacC%eps*(1._dp + rN%BacR%aA1))*rN%BacR%r(1) &
              + (1._dp + rN%SRBC%eps*(rN%SRBR%n1 - 1._dp))*rN%SRBR%r(1) + rN%PhC%eps*(1._dp - rN%PhR%aA1)*rN%PhR%r(1) &
              - rN%SOxC%eps*rN%SOxR%r(1)
      ! inorganic phosphate, h3po4
      fval(4) = transport(4) + BioS_P*( -rN%PhyC%eps*rN%PhyR%r(2) - rN%GSBC%eps*rN%GSBR%r(2) + (1._dp - rN%GzC%eps)**2*sum(rN%GzR%r(1:8)) &
              + (1._dp - rN%AGzC%eps)**2*sum(rN%AGzR%r(1:8)) - rN%BacC%eps*rN%BacR%r(1) - rN%SRBC%eps*rN%SRBR%r(1) - rN%PhC%eps*rN%PhR%r(1) ) + rN%BacR%r(3) + rN%SRBR%r(3) &
              - rN%SOxC%eps*BioS_P*rN%SOxR%r(1)
      ! sulfate, h2so4
      fval(5) = transport(5) + 0.5_dp*rN%GSBC%eps*rN%GSBR%r(1) - (0.5_dp + rN%GSBC%eps*(rN%GSBR%aA2 + 0.5_dp*(rN%GSBR%n(2) - 1._dp)))*rN%GSBR%r(2) &
              - rN%AGzR%aCi*(1._dp - rN%AGzC%eps)*sum(rN%AGzR%r(1:8)) - (0.5_dp + rN%SRBC%eps*(rN%SRBR%aA1 + 0.5_dp*(rN%SRBR%n1 - 1._dp)))*rN%SRBR%r(1) &
              + (rN%SOxC%eps*(rN%SOxR%aA1 + rN%SOxR%n1 - 1._dp) + 1._dp)*rN%SOxR%r(1)
      ! Hydrogen sulfide, h2s
      fval(6) = transport(6) - 0.5_dp*rN%GSBC%eps*rN%GSBR%r(1) + (0.5_dp + rN%GSBC%eps*(rN%GSBR%aA2 + 0.5_dp*(rN%GSBR%n(2) - 1._dp)))*rN%GSBR%r(2) &
              + rN%AGzR%aCi*(1._dp - rN%AGzC%eps)*sum(rN%AGzR%r(1:8)) + (0.5_dp + rN%SRBC%eps*(rN%SRBR%aA1 + 0.5_dp*(rN%SRBR%n1 - 1._dp)))*rN%SRBR%r(1) &
              - (rN%SOxC%eps*(rN%SOxR%aA1 + rN%SOxR%n1 - 1._dp) + 1._dp)*rN%SOxR%r(1)
      ! Carbon stored by Phy, C_Phy
      fval(7) = transport(7) + rN%PhyC%eps*rN%PhyR%r(1) - (1._dp + rN%PhyC%eps*rN%PhyR%n(2))*rN%PhyR%r(2) - Cn%C_Phy/Cn%PhyS*(rN%GzR%r(1) + rN%AGzR%r(1))     
      ! Carbon stored by GSB, C_GSB
      fval(8) = transport(8) + rN%GSBC%eps*rN%GSBR%r(1) - (1._dp + rN%GSBC%eps*rN%GSBR%n(2))*rN%GSBR%r(2) - Cn%C_GSB/Cn%GSBS*(rN%GzR%r(2) + rN%AGzR%r(2))
      ! Labile carbon, C_L
      fval(9) = transport(9) + Cn%C_Phy/Cn%PhyS*(rN%GzR%r(1) + rN%AGzR%r(1)) + Cn%C_GSB/Cn%GSBS*(rN%GzR%r(2) + rN%AGzR%r(2)) - rN%BacR%r(1) + rN%BacR%r(2) &
              - (1._dp + rN%SRBC%eps*rN%SRBR%n1)*rN%SRBR%r(1) + rN%SRBR%r(2) - rN%PhC%eps*rN%PhR%r(1)
      ! Detrital carbon, C_D
      fval(10) = transport(10) + rN%GzC%eps*(1._dp - rN%GzC%eps)*sum(rN%GzR%r(1:8)) + rN%AGzC%eps*(1._dp - rN%AGzC%eps)*sum(rN%AGzR%r(1:8)) - rN%BacR%r(2) - rN%SRBR%r(2)
      ! Deitrital phosphate, P_D
      fval(11) = transport(11) + BioS_P*( rN%GzC%eps*(1._dp - rN%GzC%eps)*sum(rN%GzR%r(1:8)) + rN%AGzC%eps*(1._dp - rN%AGzC%eps)*sum(rN%AGzR%r(1:8)) ) - rN%BacR%r(3) - rN%SRBR%r(3)
      ! Phy biological structure, PhyS
      fval(12) = transport(12) + rN%PhyC%eps*rN%PhyR%r(2) - rN%GzR%r(1) - rN%AGzR%r(1)
      ! GSB biological structure, GSBS
      fval(13) = transport(13) + rN%GSBC%eps*rN%GSBR%r(2) - rN%GzR%r(2) - rN%AGzR%r(2)
      ! Grazer biological structure, GzS
      fval(14) = transport(14) + rN%GzC%eps*sum(rN%GzR%r(1:8)) - rN%GzR%r(3) - rN%AGzR%r(3)
      ! Anaerobic grazer biological structure, AGzS
      fval(15) = transport(15) + rN%AGzC%eps*sum(rN%AGzR%r(1:8)) - rN%GzR%r(4) - rN%AGzR%r(4)
      ! Bacterial biological structure, BacS
      fval(16) = transport(16) + rN%BacC%eps*rN%BacR%aA1*rN%BacR%r(1) - rN%GzR%r(5) - rN%AGzR%r(5)
      ! Sulfate reducing bacteria biological structure, SRBS
      fval(17) = transport(17) + rN%SRBC%eps*rN%SRBR%r(1) - rN%GzR%r(6) - rN%AGzR%r(6)
      ! Photoheterotroph biological structure, PhS
      fval(18) = transport(18) + rN%PhC%eps*rN%PhR%aA1*rN%PhR%r(1) - rN%GzR%r(7) - rN%AGzR%r(7)
      ! Sulfur bacteria, SOxS
      fval(19) = transport(19) + rN%SOxC%eps*rN%SOxR%r(1) - rN%GzR%r(8) - rN%AGzR%r(8)

      ! Now, remove the dimensionality of the state vector, except for time.
      fval(1:npde) = fval(1:npde)/Cc(1:npde)      
      ! Return calling values to their initial dimensionless forms
      !x = x/Xc
      !c(1:npde) = c(1:npde)/Cc(1:npde)
      !cx(1:npde) = cx(1:npde)/Cc(1:npde)*Xc
      !cxx(1:npde) = cxx(1:npde)/Cc(1:npde)*Xc*Xc      
      return           
   end subroutine PDEs

   subroutine bndxa(t, u, ux, bval, npde)
      ! Boundary conditions at left or surface of pond.
      use realPrec
      use parameters
      use typeVars      
      ! Thermo module gives solubility of O2 (microMolar/atm), CO2 (uM/atm) and H2S in water and requires temp (K), ionic strenght (M), and pH
      ! For CO2, the Henry's value returns the total DIC that would be in equilibrium with the CO2 at a given atm.
      use thermoData, only: solO2, co2g2l, freeCO2, solH2S 
      use functions, only: stepUp
      implicit none
      !-----------------------------------------------------------------------
      ! PURPOSE:
      !       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
      !       LEFT SPATIAL END POINT X = XA.
      !                           B(T, c, cx) = 0
      !
      !-----------------------------------------------------------------------
      ! SUBROUTINE PARAMETERS:
      ! INPUT:
      real(dp), intent(in)   ::  t           ! time (d)
      real(dp), intent(inout)::  u(npde)     ! value of c at (t,x)
      real(dp), intent(inout)::  ux(npde)    ! first spatial derivative (dc/dx)
      real(dp), intent(out)  ::  bval(npde)  ! Boundary condition at the left boundary point (i.e., x = 0, at surface)
      integer,  intent(in)   ::  npde        ! number of PDEs
      ! local declarations
      type(concentrations) Cn ! used to give names to state vars.      
      real(dp) D, dDdx, ccA, dadx, T_K, pH, is, o2Eq, co2Eq, co2, h2sEq
      real(dp) c(maxPDE), cx(maxPDE)
      real(dp) vS(maxPDE)
      real(dp) spinUp

      ! Convert dimensionales quantities to their dimensional forms
      ! All calculations will be in demensional quantities, then non dimensionalized on return.
      c = u*Cc(1:npde)
      cx = ux*Cc(1:npde)/Xc      
      where (c < 0._dp) c = absZero ! prevent negative terms
      
      ! Copy state variables values to useful names
      call setCnames (c(1:nconc), c(nconc+1:nconc+nbioS), Cn)   
      ! Get or calculate needed variables
      call dispersion (t, 0._dp, D, dDdx) ! Dispersion and it's derivative (m2/d) at surface
      call ccArea     (t, 0._dp, ccA, dadx) ! crossectional area (m2) at surface     
      call tempK (t, 0._dp, T_K) ! get temperature at surface at time t
      call pHtx (t, 0._dp, pH) ! pH at surface at time t
      spinUp = stepUp(t, tStep, tSig) ! This allows easier starting
      vS = vSink*spinUp

      is = 0.72*Cn%sal/35._dp  ! the ionic strength of seawater is 0.72 M, so estimate ionic strenght by ratio of salinities.
      
      ! The flux at any location in the PDE domain, including the boundaries, is given by
      ! - D*cx(1:npde) + q*c(1:npde)/ccA + vS(1:npde)*c(1:npde)
      ! This equation is then set to what ever that specified boundary flux is.
      ! For the surface BC, set the flux equal to q*c(1:npde)/ccA, that is
      ! q*c(1:npde)/ccA = - D*cx(1:npde) + q*c(1:npde)/ccA + vS(1:npde)*c(1:npde)
      ! Or this can be rewritten as
      ! - D*cx(1:npde) + vS(1:npde)*c(1:npde) = 0
      ! So, B can be set to the LHS of the above equation
      bval(1:npde) = - D*cx(1:npde) + vS(1:npde)*c(1:npde) ! This is the dimensional version, but divided by ccA
      ! Only the O2, DIC and H2S state variables need to be modified to account for air-water gas exchange            
      ! ** Flux into water is positive **
      ! For oxygen
      o2Eq = pO2*solO2(T_K, is, pH) ! O2 in equilbrium with atmosphere O2 pressure set in pO2
      bval(2) = bval(2) - pV_o2 *( o2Eq - Cn%o2 )*spinUp ! O2 flux
      ! for CO2
      co2Eq = pCO2*co2g2l(T_K, is, pH)
      co2 = Cn%dic*freeCO2(T_K, is, pH) ! Calculate free CO2 in water from DIC
      bval(3) = bval(3) - pV_co2*( co2Eq - co2 )*spinUp !CO2 flux 
      ! for H2S
      h2sEq = pH2S*solH2S(T_K, is, pH) ! H2S in equilbrium with atmosphere H2S pressure (atm) set in pH2S
      bval(6) = bval(6) - pV_h2s *( h2sEq - Cn%h2s )*spinUp ! H2S flux
      
      ! Now, remove the dimensionality of the boundary vector.  
      ! Note, in MEP-phototroph-1D-OC_V1.14.f90, I only removed the concentration dimensions, but the units
      ! of bval as shown above are mmol/m2/d.  Since time is not been removed from the problem, the units that should be removed are mmol/m2
      ! so bval should be divided by Cc*Xc, not just Cc. It would be good to rerun V1.14 to see if this matters much, because ultimately bval is driven to zero.
      ! V1.14 was rerun (as 1.0b) with the additional removal of length by Xc, and this did not have any impact.
      bval(1:npde) = bval(1:npde)/(Cc(1:npde)*Xc)      
      ! Return calling values to their initial dimensionless forms
      !c(1:npde) = c(1:npde)/Cc(1:npde)
      !cx(1:npde) = cx(1:npde)*Xc/Cc(1:npde)      
      
      return
   end subroutine bndxa

   subroutine bndxb(t, u, ux, bval, npde)
      ! This routine is to specify the right or bottom boundary conditions.
      use realPrec
      use parameters
      use typeVars, only: concentrations, maxPDE
      use functions, only: stepUp
      implicit none
      !-----------------------------------------------------------------------
      ! PURPOSE:
      !       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
      !       RIGHT SPATIAL END POINT X = XB.
      !                           B(T, U, UX) = 0
      !
      !-----------------------------------------------------------------------
      ! SUBROUTINE PARAMETERS:
      real(dp), intent(in)   :: t           ! time (d)
      real(dp), intent(inout):: u(npde)     ! value of c at (t,x)
      real(dp), intent(inout):: ux(npde)    ! first spatial derivative (dc/dx)
      real(dp), intent(out)  :: bval(npde)  ! Boundary condition at the right boundary point (i.e., x = depth, at bottom)
      integer , intent(in)   :: npde        ! number of PDEs
      ! local declarations
      type(concentrations) Cn ! used to give names to state vars.      
      real(dp) flux(maxPDE), D, dDdx, ccA, dadx, q, qL, dqdx, cL(maxPDE)
      real(dp) f_PC, f_PH, f_PO, f_PN, f_PP, alp_P, bet_P, gam_P, del_P, aO2, bO2, aSO4, bSO4, flux_O2, flux_SO4, flux_DIC, flux_PO4
      real(dp) c(maxPDE), cx(maxPDE)
      real(dp) vS(maxPDE)
      real(dp) spinUp

      ! Convert dimensionales quantities to their dimensional forms
      ! All calculations will be in demensional quantities, then non dimensionalized before returning.
      c = u*Cc(1:npde)
      cx = ux*Cc(1:npde)/Xc
      where (c < 0._dp) c = absZero ! prevent negative terms
            
      ! Copy state variables values to useful names
      call setCnames (c(1:nconc), c(nconc+1:nconc+nbioS), Cn)            
      ! get transport terms
      call dispersion (t, depth, D, dDdx) ! Dispersion and it's derivative (m2/d)
      call ccArea     (t, depth, ccA, dadx) ! crossectional area (m2)
      call volFlow    (t, depth, q, dqdx, qL, cL) ! volumetric flow rate (m3/d) and laterial inputs if any. qL is in m2/d
      spinUp = stepUp(t, tStep, tSig) ! This allows easier starting
      q = q*spinUp
      dqdx = spinUp*dqdx
      qL = qL*spinUp     
      vS = vSink*spinUp ! this slowly spins up sinking.
      ! flux here represents what the acutal flux will be from the model.  What it is equal to must be set to determine bval
      flux(1:npde) = - D*cx(1:npde) + q*c(1:npde)/ccA + vS(1:npde)*c(1:npde) ! This is the dimensional version, but divided by ccA
      
      ! BC at sediment interface. Here their is a flow in of some of the dissolved consituents, so the concentration
      ! in that incoming water needs to be set by concBC and bioSBC.  
      
      associate ( sal_BC => concBC(1,2), o2_BC => concBC(2,2), dic_BC => concBC(3,2), po4_BC => concBC(4,2), so4_BC => concBC(5,2), h2s_BC => concBC(6,2), &
                  C_Phy_BC => concBC(7,2), C_GSB_BC => concBC(8,2), C_L_BC => concBC(9,2), C_D_BC => concBC(10,2), P_D_BC => concBC(11,2), &
                  phyS_BC => bioSBC(1,2), GSBS_BC => bioSBC(2,2), GzS_BC => bioSBC(3,2), AGzS_BC => bioSBC(4,2), BacS_BC => bioSBC(5,2), SRBS_BC => bioSBC(6,2), &
                  PhS_BC => bioSBC(7,2), SOxS_BC => bioSBC(8,2) )       
         
         ! calculate the flux of sinking material to the sediments, as this will be used to set the BC for various state variables.  
         ! First determine the bulk C, H, O, N, P composition of the sinking material.  f_PC is the carbon flux
         f_PC = vS(7)*Cn%C_Phy + vS(8)*Cn%C_GSB + vS(10)*Cn%C_D &
              + vS(12)*Cn%PhyS + vS(13)*Cn%GSBS + vS(14)*Cn%GzS + vS(15)*Cn%AGzS + vS(16)*Cn%BacS + vS(17)*Cn%SRBS + vS(18)*Cn%PhS + vS(19)*Cn%SOxS
         ! There is only a tiny amount of H comming with P_D, so do not include with f_PH
         f_PH = 2._dp*(vS(7)*Cn%C_Phy + vS(8)*Cn%C_GSB + vS(10)*Cn%C_D) &
              + BioS_H*(vS(12)*Cn%PhyS + vS(13)*Cn%GSBS + vS(14)*Cn%GzS + vS(15)*Cn%AGzS + vS(16)*Cn%BacS + vS(17)*Cn%SRBS + vS(18)*Cn%PhS + vS(19)*Cn%SOxS)
         f_PO = vS(7)*Cn%C_Phy + vS(8)*Cn%C_GSB + vS(10)*Cn%C_D + 4._dp*vS(11)*Cn%P_D &
              + BioS_O*(vS(12)*Cn%PhyS + vS(13)*Cn%GSBS + vS(14)*Cn%GzS + vS(15)*Cn%AGzS + vS(16)*Cn%BacS + vS(17)*Cn%SRBS + vS(18)*Cn%PhS + vS(19)*Cn%SOxS)
         f_PN = BioS_N*(vS(12)*Cn%PhyS + vS(13)*Cn%GSBS + vS(14)*Cn%GzS + vS(15)*Cn%AGzS + vS(16)*Cn%BacS + vS(17)*Cn%SRBS + vS(18)*Cn%PhS + vS(19)*Cn%SOxS)
         f_PP = vS(11)*Cn%P_D + BioS_P*(vS(12)*Cn%PhyS + vS(13)*Cn%GSBS + vS(14)*Cn%GzS + vS(15)*Cn%AGzS + vS(16)*Cn%BacS + vS(17)*Cn%SRBS + vS(18)*Cn%PhS + vS(19)*Cn%SOxS)
         ! Determine C-normalized composition of sinking POM
         alp_P = f_PH/f_PC
         bet_P = f_PO/f_PC
         gam_P = f_PN/f_PC
         del_P = f_PP/f_PC
         ! The aerobic decompoistion of sinking POM is
         ! C Halp_p OBet_P Ngam_P Pdel_P + aO2 -> H2CO3 + gam_P NH3 + del_P H3PO4 + bO2 H2O
         aO2 = (4._dp + alp_P - 2._dp*bet_P - 3._dp*gam_P + 5._dp*del_P)/4._dp
         bO2 = (-2._dp + alp_P - 3._dp*gam_P - 3._dp*del_P)/4._dp
         ! The anaerobic decompoistion of sinking POM with H2SO4 is
         ! C Halp_p OBet_P Ngam_P Pdel_P + aSO4 H2SO4 -> H2CO3 + gam_P NH3 + del_P H3PO4 + bSO4 H2O + aSO4 H2S
         aSO4 = (4._dp + alp_P - 2._dp*bet_P - 3._dp*gam_P + 5._dp*del_P)/8._dp
         bSO4 = (-2._dp + alp_P - 3._dp*gam_P - 3._dp*del_P)/4._dp
         ! The fluxes of O2, DIC, H2SO4 and H2S due to sinking POM decomposition are given by the following, where it is assumed that first O2 is used, followed by
         ! H2SO4, and any POC remaining is fermented to CO2 + CH4.  
         ! Note, the flux of material to the benthose can cause problems in initializing the PDE solver (i.e., problem getting consistent BC and IC)
         ! To avoid this problem, the flux to the benthose can be steped up gracefully with the stepUp function
         ! f_PC = f_PC*stepUp(t, tStep, tSig) ! tStep and tsig are parameters in the *.inp file.  THIS IS NOW HANDLED VIA MODIFICATION TO vS Instead.
         flux_O2  = aO2*f_OMremin*f_PC*Cn%o2/(Cn%o2 + k_O2remin) ! this is the oxygen that is consumed by the particulate flux.  It may not consume all of f_OMremin*f_PC
         flux_SO4 = aSO4*(f_OMremin*f_PC - flux_O2/aO2)*Cn%so4/(Cn%so4 + k_SO4remin) ! SO4 only is used of there is not enough O2 to consume all the carbon flux
         ! Whatever sinking carbon is fermented to CH4 + CO2, given by 1/2(f_OMremin*f_PC - flux_O2/aO2 - flux_SO4/aSO4), but only tracking DIC, not CH4, so
         flux_DIC = 0.5_dp*(f_OMremin*f_PC + flux_O2/aO2 + flux_SO4/aSO4) ! this accounts for CO2 coming from methanogenesis also.
         flux_PO4 = del_P*f_OMremin*f_PC ! all the P in the sining PC is remineralized.                 
         ! salt.  
         !Here, flux is set by the incoming water (q is < 0) and the salt in it.
         bval(1) = flux(1) - (q*sal_BC/ccA + vS(1)*Cn%sal) ! not, for dissolved consituents, vS(i) will be zero, but include here and below for kicks.
         ! O2.  
         ! In addition to the O2 in the incoming water, there is also a sediment O2 demand out of the water column.  
         !  This is determined based on the OM flux that is leaving the water column by sinking and can be oxidized.
         bval(2) = flux(2) - (q*o2_BC/ccA + vS(2)*Cn%o2 + flux_O2)
         ! DIC
         ! In addition to the DIC transported in by q, there is also SOURCE of DIC retruned by OM degridation
         bval(3) = flux(3) - (q*dic_BC/ccA + vS(3)*Cn%dic - flux_DIC) 
         ! PO4
         ! This is similar to DIC
         bval(4) = flux(4) - (q*po4_BC/ccA + vS(4)*Cn%po4 - flux_PO4)          
         ! so4
         bval(5) = flux(5) - (q*so4_BC/ccA + vS(5)*Cn%so4 + flux_SO4)          
         ! h2s
         bval(6) = flux(6) - (q*h2s_BC/ccA + vS(6)*Cn%h2s - flux_SO4)          
         ! C_Phy
         ! Here, flux is set by what is coming in with water (will be set to zero usually)
         bval(7) = flux(7) - (q*C_Phy_BC/ccA + vS(7)*Cn%C_Phy)
         ! C_GSB
         bval(8) = flux(8) - (q*C_GSB_BC/ccA + vS(8)*Cn%C_GSB)
         ! C_C_L
         bval(9) = flux(9) - (q*C_L_BC/ccA + vS(9)*Cn%C_L)
         ! C_D (detritus)
         bval(10) = flux(10) - (q*C_D_BC/ccA + vS(10)*Cn%C_D)
         ! P_D (detritus)
         bval(11) = flux(11) - (q*P_D_BC/ccA + vS(11)*Cn%P_D)
         ! bS_a
         ! All biological structures are similar.
         ! PhyS
         bval(12) = flux(12) - (q*PhyS_BC/ccA + vS(12)*Cn%PhyS)
         ! GSBS
         bval(13) = flux(13) - (q*GSBS_BC/ccA + vS(13)*Cn%GSBS)
         ! GzS
         bval(14) = flux(14) - (q*GzS_BC/ccA + vS(14)*Cn%GzS)
         ! AGzS
         bval(15) = flux(15) - (q*AGzS_BC/ccA + vS(15)*Cn%AGzS)
         ! BacS
         bval(16) = flux(16) - (q*BacS_BC/ccA + vS(16)*Cn%BacS)
         ! SRBS
         bval(17) = flux(17) - (q*SRBS_BC/ccA + vS(17)*Cn%SRBS)
         ! PhS
         bval(18) = flux(18) - (q*PhS_BC/ccA + vS(18)*Cn%PhS)
         ! SOxS
         bval(19) = flux(19) - (q*SOxS_BC/ccA + vS(19)*Cn%SOxS)
      end associate
      
      ! Now, remove the dimensionality of the boundary vector.  
      ! Note, in MEP-phototroph-1D-OC_V1.14.f90, I only removed the concentration dimensions, but the units
      ! of bval as shown above are mmol/m2/d.  Since time is not been removed from the problem, the units that should be removed are mmol/m2
      ! so bval should be divided by Cc*Xc, not just Cc. It would be good to rerun V1.14 to see if this matters much, because ultimately bval is driven to zero.
      bval(1:npde) = bval(1:npde)/(Cc(1:npde)*Xc)
      ! Return calling values to their initial dimensionless forms
      !c(1:npde) = c(1:npde)/Cc(1:npde)
      !cx(1:npde) = cx(1:npde)*Xc/Cc(1:npde)      
      
      return
   end subroutine bndxb

   subroutine Cinit(X, C, NPDE)
      use realPrec
      use parameters
      implicit none
      !-----------------------------------------------------------------------
      ! PURPOSE:
      !       THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL 
      !       CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0 
      !       AT THE SPATIAL COORDINATE X.
      !-----------------------------------------------------------------------
      ! SUBROUTINE PARAMETERS:
      ! INPUT:
      DOUBLE PRECISION        X
      !                               THE SPATIAL COORDINATE (dimensionless).
      !
      INTEGER                 NPDE
      !                               THE NUMBER OF PDES IN THE SYSTEM.
      !
      ! OUTPUT:
      DOUBLE PRECISION        C(NPDE)
      !                               C(1:NPDE) IS VECTOR OF INITIAL VALUES OF
      !                               THE UNKNOWN FUNCTION AT T = T0 AND THE
      !                               GIVEN VALUE OF X. (Dimensionless)
      !-----------------------------------------------------------------------
      ! local declarations
      real(dp) delt, xD
      integer nlow
   
      ! these are deterined from what was stored at the end of the previous interval
      ! last state variable is internal entropy produced.
      ! Previous solution was saved on Xpts grid
      ! Interpolate time in tOCnodes
      ! convert to dimensional units for x
      xD = x*Xc
      call interp1D (xD, nXpts, Xpts, nlow, delt)
      ! concIni and bioSini have full units
      C(1:nconc) = concIni(1:nconc,nlow) + delt*( concIni(1:nconc,nlow+1) - concIni(1:nconc,nlow) )
      C(nconc+1:npde) = bioSini(1:nbioS,nlow) + delt*( bioSini(1:nbioS,nlow+1) - bioSini(1:nbioS,nlow) )
      C = C/Cc(1:npde) ! nondimensionalize for BACOLI
      return
   end subroutine Cinit

   subroutine dispersion (t, x, D, dDdx)
      ! Calculate Dispersion (m2/d) and it's derivative at time t (d) and location x (m)
      use realPrec
      use parameters
      implicit none
      real(dp), intent(in)  ::  t           ! time (d)
      real(dp), intent(in)  ::  x           ! spatial coord (m)
      real(dp), intent(out) ::  D           ! Dispersion at (t,x)
      real(dp), intent(out) ::  dDdx        ! first spatial derivative (dD/dx)
      ! local declarations
      D = 0.09981082055826766 - 0.006533097434132809*x + 0.0009552253162313702*x**2 + 0.0006952972126051994*x**3 &
         - 0.00022146974772627084*x**4 + 0.000018962538858758408*x**5 - 4.960920581257112e-7*x**6
      dDdx = -0.006533097434132809 + 0.0019104506324627403*x + 0.0020858916378155984*x**2 - 0.0008858789909050833*x**3 + 0.00009481269429379204*x**4 - 2.9765523487542675e-6*x**5
      return
   end subroutine dispersion
   
   subroutine ccArea (t, x, a, dadx)
      ! Calculate cross sectional area and its first spatial derivative at (t, x)
      use realPrec
      use parameters
      implicit none
      real(dp), intent(in)  ::  t           ! time (d)
      real(dp), intent(in)  ::  x           ! spatial coord (m)
      real(dp), intent(out) ::  a           ! area (m2) at (t,x)
      real(dp), intent(out) ::  dadx        ! first spatial derivative (da/dx)
      ! local declarations
      a = 146372.78791478608_dp - 23695.92113837067_dp*x + 3151.809111431212_dp*x**2 - 285.1969685799414_dp*x**3 + 9.159821439465391_dp*x**4
      dadx =-23695.92113837067_dp + 6303.618222862424_dp*x - 855.5909057398243_dp*x**2 + 36.639285757861565_dp*x**3
      return
   end subroutine ccArea
   
   subroutine volFlow (t, x, q, dqdx, qL, cL)
      ! calculates the volumetric flow rate (m3/d) and it's first spatial derivative
      ! also provided is the laterial inputs qLat (m2/d) and concentration of state variables in that flow
      use realPrec
      use parameters
      implicit none
      real(dp), intent(in) :: t              ! time (d)
      real(dp), intent(in) :: x              ! spatial coord (m)
      real(dp), intent(out):: q              ! volumetric flow (m3/d) at (t,x)
      real(dp), intent(out):: dqdx           ! first spatial derivative (dq/dx)
      real(dp), intent(out):: qL             ! laterial input (m2/d) at (t,x)
      real(dp), intent(out):: cL(nconc+nbioS)   ! concentration of state variables in qLat, does NOT include sigma variable
      ! local declarations
      real(dp), parameter:: qBottom = -721._dp ! this is the flow entering the "bottom" that is a mix of all SW plus some FW entrainment
      real(dp), parameter:: a = 0.5_dp ! exponent (1/m)
      real(dp), parameter:: qObs = 2429._dp ! total laterial inupts of FW (m3/d), not including that entrained
      real(dp) qL0 ! laterial flow at x = 0
      real(dp) IqL ! Integrated laterial flow needed to get q
      
      qL0 = qObs*a*exp(a*depth)/(exp(a*depth) - 1._dp)
      qL = qL0*exp(-a*x)
      IqL = qL0*( exp(-depth*a) - exp(-a*x) )/a
      q = qBottom + IqL  ! qL contributes to total upward flow (in fact, it's most of it), which is all < 0
      dqdx = qL0*exp(-a*x) ! this is the derivative
      cL   = cLat  ! set by input file.
      return
   end subroutine volFlow
   
   subroutine tempK (t, x, T_K)
      ! This routine returns the temperature (K) at time t (d) and location x (m)
      use realPrec
      use parameters
      implicit none
      real(dp), intent(in) :: t    ! time (d)
      real(dp), intent(in) :: x    ! spatial coord (m)
      real(dp), intent(out):: T_K  ! temperature (K) at (t,x)
      ! These values are approximate averages over the sampling period.  Not much spread, except for +- 0.5 d
      real(dp), dimension(10), parameter:: depths = [0., 0.5, 2., 3., 4., 6., 8., 10., 12., 15.]
      real(dp), dimension(10), parameter:: T_Cobs = [24.6, 24.6, 24.4, 22.5, 18.7, 15.1, 12.6, 12.1, 12.5, 12.5] ! ends are just extrapolated flat
      real(dp) delt
      integer nlow
      call interp1D (x, 10, depths, nlow, delt)
      T_K = T_Cobs(nlow) + delt*(T_Cobs(nlow+1) - T_Cobs(nlow))
      T_K = T_K + 273.15_dp ! convert to K
      !T_K = 298._dp ! for now just set as constant
      return
   end subroutine tempK
   
   subroutine pHtx (t, x, pH)
      ! This routine returns the pH at time t (d) and location x (m)
      use realPrec
      use parameters
      implicit none
      real(dp), intent(in) :: t    ! time (d)
      real(dp), intent(in) :: x    ! spatial coord (m)
      real(dp), intent(out):: pH  ! temperature (K) at (t,x)
      ! local declarations
      ! These values are the average values for each depth (plus extrapolation) obtained from the Siders Pond 24 hr sampling.
      real(dp), dimension(10), parameter:: depths = [0., 0.5, 2., 3., 4., 6., 8., 10., 12., 15.]
      real(dp), dimension(10), parameter:: pHobs = [7.57143, 7.575714286, 7.588571429, 8.14, 7.055714286, 6.517142857, 6.547142857, 6.898571429, 6.921428571, 6.95571]
      real(dp) delt
      integer nlow
   
      call interp1D (x, 10, depths, nlow, delt)
      pH = pHobs(nlow) + delt*(pHobs(nlow+1) - pHobs(nlow))
      !pH = 7.6_dp ! for now just set to 7
      return
   end subroutine pHtx   
   
      
   
   