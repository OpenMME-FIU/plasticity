#include "../../../include/crystalPlasticity.h"
#include "userFunctions.cc"
#include "CP_UMAT_FCC.cpp"
//////////////////////////////////////////////////////////////////////////
//////calculatePlasticity.cc numerically integrates the constitive model.
//////This calculatePlasticity.cc is based on the following rate-dependent crystal plasticity model:
//////SR Kalidindi, Polycrystal plasticity: constitutive modeling and deformation processing,
////// PhD thesis, MIT, 1992. In addition to isotropic hardening tehre is also kinematic hardening
////// where the backstress evolves based on the OW hardening model. The guess stress to start the
////// Newton-Raphson iteration is the previously converged stress. Exponential update of Fp is
//// implemented along with a line-search algorithm to solve the nonlinear system. The
//// tangent modulus computation is also more involved.
//////To use this file, one should copy this calculatePlasticity.cc into the following folder:
//////    plasticity/src/materialModels/crystalPlasticity/
////// (It should replace the original calculatePlasticity.cc inside that folder)
//////Next, they should copy the file userFunctions.cc into the following folder:
//////    plasticity/src/materialModels/crystalPlasticity/
//////
//////Finally, the PRISMS-Plasticity should be recompiled.
//////////////////////////////////////////////////////////////////////////////
////
//
template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
  unsigned int quadPtID,
  unsigned int StiffnessCalFlag)
  {

    multiphaseInit(cellID,quadPtID);



    /////////////////////////////////////////////////////////
    FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim),F_t(dim,dim);  //Elastic, Plastic, and total deformation gradient
    FullMatrix<double> E_t(dim,dim) // Green Strain for UMAT
    Vector<double> s_alpha_t(n_Tslip_systems),slipfraction_t(n_slip_systems),twinfraction_t(n_twin_systems); // Slip resistance
    Vector<double> W_kh_t(n_Tslip_systems),W_kh_t1(n_Tslip_systems),W_kh_t2(n_Tslip_systems),signed_slip_t(n_Tslip_systems) ;
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)
    FullMatrix<double> Tinter_diff_guess(dim,dim), Ep_t (dim,dim), iMinusL (dim,dim);
    double Ep_eff_cum_t;
    double tol1, delgam_ref,strexp,locres_tol,locres_tol2,h_back1,h_back2,r_back1,r_back2,m_back1,m_back2,back_lim_1,back_lim_2,b_back1,b_back2;
    unsigned int nitr1,nitr2;
    unsigned int ii;
    FullMatrix<double> rotmat(dim,dim); // Rotation matrix of the crystal orientation
    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim),temp4(dim,dim) ; // Temporary matrices
    FullMatrix<double> temp5(dim,dim),temp6(dim,dim), temp7(dim,dim), temp8(dim,dim); // Temporary matrices
    FullMatrix<double> T_tau(dim,dim),P_tau(dim,dim); // Stress measures for current timestep
    FullMatrix<double> FP_inv_t(dim,dim),FE_tau_trial(dim,dim),CE_tau_trial(dim,dim) ; // Kinematic descriptors
    FullMatrix<double>FP_inv_tau(dim,dim),F_inv_tau(dim,dim); // Kinematic descriptors
    FullMatrix<double> SCHMID_TENSOR1(n_Tslip_systems*dim,dim),Normal_SCHMID_TENSOR1(n_Tslip_systems*dim,dim); // Projection matrices
    Vector<double> m1(dim),n1(dim);

    // Elastic Modulus
    FullMatrix<double> Dmat(2*dim,2*dim),Dmat2(2*dim,2*dim),TM(dim*dim,dim*dim);
    Vector<double> vec2(dim*dim);
    Vector<double> s_alpha_tau(n_Tslip_systems),slipfraction_tau(n_slip_systems),twinfraction_tau(n_twin_systems) ;
    Vector<double> W_kh_tau(n_Tslip_systems),W_kh_tau1(n_Tslip_systems),W_kh_tau2(n_Tslip_systems); // Converged backstresses
    Vector<double> W_kh_tau_it(n_Tslip_systems),W_kh_tau1_it(n_Tslip_systems),W_kh_tau2_it(n_Tslip_systems); // Iterative backstresses
    Vector<double> h_beta(n_Tslip_systems),h0(n_Tslip_systems),a_pow(n_Tslip_systems),s_s(n_Tslip_systems); // Isotropic hardening parameters
    FullMatrix<double> h_alpha_beta_t(n_Tslip_systems,n_Tslip_systems),Ep_tau(dim,dim),del_Ep_tau(dim,dim);
    Vector<double> resolved_shear_tau(n_Tslip_systems), normal_stress_tau(n_Tslip_systems),signed_slip_tau(n_Tslip_systems);
    double Ep_eff_cum_tau,del_Ep_eff_cum_tau;
    FullMatrix<double> PK1_Stiff(dim*dim,dim*dim); // Tangent modulus
    FullMatrix<double> T_star_tau(dim,dim),mtemp(dim,dim);
    Vector<double> vtemp(2*dim);
    double det_FE_tau, det_F_tau, trny_op, trny_in;
    double m1_norm, n1_norm ;

    unsigned int itr1, itr2;
    double sctmp1, sctmp2, sctmp3, sctmp4; // Scalars used as temporary variables
    Vector<double> G_iter(2*dim),locres_vec(2*dim+2*n_Tslip_systems),stateVar_it(2*dim+2*n_Tslip_systems),stateVar_temp(2*dim+2*n_Tslip_systems),stateVar_diff(2*dim+2*n_Tslip_systems),T_star_iter_vec(2*dim);
    Vector<double> gradFold(2*dim+2*n_Tslip_systems) ;
    FullMatrix<double> btemp1(2*dim,2*dim),J_iter(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems),J_iter_inv(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems),LP_acc(dim,dim),LP_acc2(dim,dim) ;
    FullMatrix<double> J_iter_cp(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems);
    FullMatrix<double> T_star_iter(dim,dim),T_star_iterp(dim,dim);
    FullMatrix<double> CE_t(dim,dim);
    Vector<double> s_alpha_it(n_Tslip_systems),s_alpha_iterp(n_Tslip_systems),delgam_tau(n_Tslip_systems),delgam_tau_iter(n_Tslip_systems);

    Vector<double> vtmp1(2*dim),vtmp2(2*dim),vtmp3(2*dim);
    Vector<double> vtmp4(n_Tslip_systems);
    Vector<double> nv1(2*dim),nv2(2*dim);
    double locres, locres2;
    double sgnm, sgnm2 , sgnm3 , sgnm4, Fold ;

    FullMatrix<double>  cnt1(dim*dim,dim*dim),cnt2(dim*dim,dim*dim),cnt3(dim*dim,dim*dim),cnt4(dim*dim,dim*dim); // Variables to track individual contributions
    FullMatrix<double> dFedF(dim*dim,dim*dim),dFpdFe(dim*dim,dim*dim); // Meaningful variables
    FullMatrix<double>  ntemp1(dim*dim,dim*dim),ntemp2(dim*dim,dim*dim),ntemp3(dim*dim,dim*dim),ntemp4(dim*dim,dim*dim),ntemp5(dim*dim,dim*dim),ntemp6(dim*dim,dim*dim); // Temporary variables
    FullMatrix<double> Pmat(n_Tslip_systems,n_Tslip_systems),Qmat(n_Tslip_systems,n_Tslip_systems),Pmat_inv(n_Tslip_systems,n_Tslip_systems),Rmat(n_Tslip_systems,n_Tslip_systems) ;
    FullMatrix<double> Smat(n_Tslip_systems,n_Tslip_systems) ;
    FullMatrix<double> Qmat2(n_Tslip_systems,dim*dim),dgamdFe(n_Tslip_systems,dim*dim) ;
    Vector<double> nvec1(dim*dim),nvec2(dim*dim) ;

    double modifier1_num, modifier1_den, modifier2_num, modifier2_den, modifier;

    double mulfac;
    FullMatrix<double> L(dim, dim);

    std::vector<double> local_twin;
    std::vector<double>::iterator result;
    double twin_pos, twin_max;
    Vector<double> quat1(4), rod(3), quat2(4), quatprod(4);

    double Criteria_Delta_F,inverseNumberOfCuts,Criteria_Delta_F_2;
    unsigned int  numberOfCuts;
    FullMatrix<double> Delta_F,div_Delta_F;
      this->pcout << "calculatePlasticity\n";

      fem::real_star_8 sse;
      fem::real_star_8 const spd = 0.0f; //
      fem::real_star_8 const scd = 0.0f;//
      fem::real_star_8 const rpl= 0.0f;//

      fem::real_star_8 const dtime= 0.0f;
      fem::real_star_8 const temp11= 0.0f;//
      fem::real_star_8 const dtemp= 0.0f;//
      fem::arr_1d<1, fem::real_star_8> predef;//
      fem::arr_1d<1, fem::real_star_8> dpred;//
      std::string cmname;//
      int const& ndi=dim;
      int const& nshr=(dim*dim - dim)/2;
      int const& ntens=dim + (dim*dim - dim)/2;
      int const& nstatv = this->userInputs.numberofUserMatStateVar1;

      int const& nprops = this->userInputs.numberofUserMatConstants1;
      fem::arr_1d<dim, fem::real_star_8> coords ;//
      fem::arr_2d<dim, dim, fem::real_star_8> drot ;//
      fem::real_star_8 pnewdt;
      fem::real_star_8 const celent = 0.0f;//
      fem::arr_1d<ntens, fem::real_star_8> ddsddt;//
      fem::arr_1d<ntens, fem::real_star_8> drplde;//
      fem::real_star_8 const drpldt = 0.0f;//
      fem::arr_1d<ntens, fem::real_star_8> stran ;//
      fem::arr_1d<ntens, fem::real_star_8> dstran;//
      fem::arr_1d<2, fem::real_star_8> time;
      int const noel = cellID;//
      int const npt = quadPtID;//
      fem::real_star_8 const layer = 0.0f;//
      int const kspt = 0;//
      int const kstep = 0;
      int const kinc = 0;//
      this->pcout << "calculatePlasticity2\n";
      fem::arr_1d<nprops, fem::real_star_8> props;
      fem::arr_1d<nstatv, fem::real_star_8> statev;
      fem::arr_1d<ntens, fem::real_star_8> stress; //Voigt notation
      fem::arr_2d<ntens, ntens, fem::real_star_8> ddsdde;
      fem::arr_2d<dim,dim, fem::real_star_8> dfgrd0;
      fem::arr_2d<dim,dim, fem::real_star_8> dfgrd1;
//      time(fem::dimension(2));
      for(int i=0; i<nstatv;i++){
          statev(i+1) = stateVar_conv[cellID][quadPtID][i];
      }
      for(int i=0; i<nprops;i++){
          props(i+1) = UserMatConstants[i];
      }
      for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
              dfgrd0(i+1,j+1)=F[i][j];
              sse += F[i][j]*P[i][j]/2; //missing divide by element volume???
          }
      }
      for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
              if(i == j)
                  stress(i+1) = T[i][j]; //add stran and E here?
              else
                  stress(7-i-j) = T[i][j];
          }
      }
//dfgrd1(fem::dimension(dim, dim));
    //////////////////////////////////////////////////////////
    cp_fcc::umat(stress,statev,ddsdde,sse,0,0,0,
                 0,0,0,0,0,time,dtime,0,0,0,0,"",
                 ndi,nshr,ntens,nstatv,props,nprops,
                 0,0,pnewdt,0, dfgrd0,dfgrd1,
                 0,0,0,0,kstep,0);
      for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
              F_tau[i][j]=dfgrd1(i+1,j+1);
          }
      }
      for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
              //check voigt ordering of umat vs this script
              unsigned int ib=0;
              unsigned int ia=0;
              if (i == j)
                  ia = i;
              else
                  ia = 6 - i - j;
              for(unsigned int k=0 ; k<dim ; k++){
                  for(unsigned int l=0 ; l<dim ; l++){

                      if (k == l)
                          ib = k;
                      else
                          ib = 6 - k - l;

                      dP_dF[i][j][k][l]=ddsdde(ia+1,ib+1);
                  }
              }
          }
      }
      for (int i; i<nstatv;i++){
          stateVar_iter[cellID][quadPtID][i]=statev(i+1);
      }
      for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
              if(i == j)
                  T_tau[i][j] = stress(i+1);
              else
                  T_tau[i][j] = stress(7-i-j);
          }
      }


    std::cout.precision(16);


    tol1=this->userInputs.modelStressTolerance;

      temp.reinit(dim, dim);
//      det_FE_tau = FE_tau.determinant();
//      T_star_tau.equ(1.0,T_star_iter);
//      FE_tau.mmult(temp, T_star_tau);
//      temp.equ(1.0 / det_FE_tau, temp);
//      temp.mTmult(T_tau, FE_tau);

      det_F_tau = F_tau.determinant();
      temp.invert(F_tau);
      F_inv_tau.equ(1.0,temp);
      P_tau.mTmult(T_tau, temp); //CHECK ME!! for transpose
      P_tau.equ(det_F_tau, P_tau);

//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    ///////////////Reading the previous converged state Variables////////////////
//    /////////////////////////////////////////////////////////////////////////////
//
//    ////////////////////////Rotating the Elsatic modulus from crystal to sample//////////////////
//    rotmat=0.0;
//    odfpoint(rotmat,rot1);
//    elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);
//    /////////////////////////////////////////////////////////////////////////////////////////////
//
//    ///////////////////////////Calculation of elastic stiffness matrix in different shapes required in the implementation//////////////
//    //Elastic Stiffness Matrix Dmat
//    Dmat.reinit(6,6) ;
//    Dmat = 0.0;
//    for (unsigned int i = 0;i<6;i++) {
//      for (unsigned int j = 0;j<6;j++) {
//        Dmat[i][j] = Dmat2[i][j];
//      }
//    }
//    for (unsigned int i = 0;i<6;i++) {
//      for (unsigned int j = 3;j<6;j++) {
//        Dmat[i][j] = 2 * Dmat[i][j];
//      }
//    }
//
//    vec2(0)=0;vec2(1)=5;vec2(2)=4;vec2(3)=5;vec2(4)=1;vec2(5)=3;vec2(6)=4;vec2(7)=3;vec2(8)=2;
//    for(unsigned int i=0;i<9;i++){
//      for(unsigned int j=0;j<9;j++){
//        TM[i][j]=Dmat2(vec2(i),vec2(j));
//      }
//    }
//    //////////////////////////////////////////////////////////////////


    P.reinit(dim, dim);T.reinit(dim, dim);S.reinit(dim, dim);
    P = P_tau;
    T = T_tau;
//    S = T_star_tau;
//
//    T_inter.reinit(dim,dim) ;
//    T_inter = T_star_tau ;




//    for(unsigned int i=0 ; i<n_twin_systems ; i++)
//    twinfraction_iter[cellID][quadPtID][i]=twinfraction_tau(i);
//
//    for(unsigned int i=0 ; i<n_slip_systems ; i++)
//    slipfraction_iter[cellID][quadPtID][i]=slipfraction_tau(i);

    if (this->userInputs.flagTaylorModel){
      F=F_tau; // Updating Deformation Gradient if it is Taylor model
    }
    /////// REORIENTATION Due to TWINNING ////////////////////

    if (enableTwinning){
      if (!this->userInputs.enableMultiphase){
        if (F_r > 0) {
          F_T = twinThresholdFraction + (twinSaturationFactor*F_e / F_r);
        }
        else {
          F_T = twinThresholdFraction;
        }
      }
      else{
        F_T = twinThresholdFraction;
      }

      //////Eq. (13) in International Journal of Plasticity 65 (2015) 61–84
      if (F_T > 1.0) {
        F_T = 1.0;
      }
      local_twin.resize(n_twin_systems,0.0);
      local_twin=twinfraction_iter[cellID][quadPtID];
      result = std::max_element(local_twin.begin(), local_twin.end());
      twin_pos= std::distance(local_twin.begin(), result);
      twin_max=local_twin[twin_pos];
      if (twin_conv[cellID][quadPtID] != 1.0) {
        if(F_r>0){
          if(twin_max > F_T){

            rod(0) = rot_conv[cellID][quadPtID][0];rod(1) = rot_conv[cellID][quadPtID][1];rod(2) = rot_conv[cellID][quadPtID][2];
            odfpoint(rotmat, rod);
            rod2quat(quat2, rod);
            quat1(0) = 0;
            quat1(1) = n_alpha[n_slip_systems + twin_pos][0];
            quat1(2) = n_alpha[n_slip_systems + twin_pos][1];
            quat1(3) = n_alpha[n_slip_systems + twin_pos][2];

            quatproduct(quatprod, quat2, quat1);


            quat2rod(quatprod, rod);

            odfpoint(rotmat, rod);

            rot_iter[cellID][quadPtID][0] = rod(0);rot_iter[cellID][quadPtID][1] = rod(1);rot_iter[cellID][quadPtID][2] = rod(2);
            rotnew_iter[cellID][quadPtID][0] = rod(0);rotnew_iter[cellID][quadPtID][1] = rod(1);rotnew_iter[cellID][quadPtID][2] = rod(2);
            twin_iter[cellID][quadPtID] = 1.0;
            for (unsigned int i = 0;i < n_twin_systems;i++) {
              s_alpha_iter[cellID][quadPtID][n_slip_systems + i] =100000;
            }
            Fe_iter[cellID][quadPtID]=IdentityMatrix(dim);
            Fp_iter[cellID][quadPtID]=IdentityMatrix(dim);
          }
        }
      }
    }

  }
  #include "../../../include/crystalPlasticity_template_instantiations.h"
