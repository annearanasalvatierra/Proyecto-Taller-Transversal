/*----------------------------------------------------------------------------

 Purpose:
   Numerical integration methods for ordinaray differential equations

   This module provides implemenation of the variable order variable 
   stepsize multistep method of Shampine & Gordon.
 
 Last modified:   2015/08/25   M. Mahooti
 
 Reference:

   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
   Freeman and Comp., San Francisco (1975).

----------------------------------------------------------------------------*/
#include "DEIntegAccel.h"
#include "sign_.h"
#include "Accel.h"
#include <cmath>

void DEIntegAccel(double t, double tout, double relerr, double abserr,int n_eqn, double y[6]){ 

    // maxnum = 500;
    double eps = 2.220446049250313e-16;
    double twou  = 2*eps;
    double fouru = 4*eps;

    typedef struct{
        int DE_INIT = 1;      // Restart integration
        int DE_DONE = 2;      // Successful step
        int DE_BADACC = 3;    // Accuracy requirement could not be achieved
        int DE_NUMSTEPS = 4;  // Permitted number of steps exceeded
        int DE_STIFF = 5;     // Stiff problem suspected
        int DE_INVPARAM = 6;  // Invalid input parameters
    } Destate;
    Destate DE_STATE;

    int State_ = DE_STATE.DE_INIT;
    bool PermitTOUT = true;         // Allow integration past tout by default
    int told = 0;

    // Powers of two (two(n)=2^n)
    double two[14]  = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};

    double gstr[14] = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468};

    double yy[n_eqn];    // Allocate vectors with proper dimension
    double wt[n_eqn];
    double p[n_eqn];
    double yp[n_eqn];

    for(int i=0; i<n_eqn; i++){
        yy[i] = 0;
        wt[i] = 0;
        p[i]  = 0;
        yp[i] = 0;
    }

    double phi[n_eqn][17];
    for(int i=0; i<n_eqn; i++){
        for(int j=0; j<17; j++){
            phi[i][j] = 0;
        }
    }

    double g[14]     = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double sig[14]   = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double rho[14]   = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double w[14]     = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double alpha[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double beta[14]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double v[14]     = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double psi_[14]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // while(true)

    // Return, if output time equals input time

    if (fabs(t - tout) < 1e-8){    // No integration
        return;
    }

    // Test for improper parameters
    double epsilon;
    if(relerr > abserr){
        epsilon = relerr;
    }else{
        epsilon = abserr;
    }

    if ( ( relerr <  0.0         ) ||             // Negative relative error bound
         ( abserr <  0.0         ) ||             // Negative absolute error bound
         ( epsilon    <= 0.0         ) ||         // Both error bounds are non-positive
         ( State_  >  DE_STATE.DE_INVPARAM ) ||   // Invalid status flag
         ( (State_ != DE_STATE.DE_INIT) &&  (t != told) ) ){
             State_ = DE_STATE.DE_INVPARAM;              // Set error code
             return;                                     // Exit
         }

    // On each call set interval of integration and counter for
    // number of steps. Adjust input error tolerances to define
    // weight vector for subroutine STEP.

    double del    = tout - t;
    double absdel = fabs(del);

    double tend   = t + 100.0*del;
    if (!PermitTOUT){
        tend = tout;
    }

    int nostep = 0;
    int kle4   = 0;
    bool stiff  = false;
    double releps = relerr/epsilon;
    double abseps = abserr/epsilon;

    bool start, OldPermit;
    double x, delsgn, h, aux;
    
    if  ( (State_==DE_STATE.DE_INIT) || (!OldPermit) || (delsgn*del<=0.0) ){
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        start  = true;
        x      = t;
        for(int i=0; i<n_eqn; i++){
            yy[i] = y[i];
        }                           
        delsgn = sign_(1.0, del);
        if(fouru*fabs(x) > fabs(tout-x)){
           aux = fouru*fabs(x);
        }else{
           aux = fabs(tout-x);
        }
        h      = sign_( aux, tout-x );
    }

    double yout[n_eqn];
    double ypout[n_eqn];
    double hi, ki, term, psijm1, gamma, eta, aux1, p5eps, round, sum, aux2, hold, hnew, temp1, temp2, temp3, temp4, temp5, temp6, tau, xold, erkm2, erkm1, erk, err, aux3, erkp1, aux4, aux5, aux6, aux7, absh, r;
    bool crash, phase1, nornd, success;
    int ifail, k, kold, kp1, kp2, km1, km2, ns, nsp1, nsp2, realns, im1, reali, nsm2, limit1, limit2, knew, i, ip1;

    while (true){   // Start step loop

          // If already past output point, interpolate solution and return
          if (fabs(x-t) >= absdel){ 
              for(int i=0; i<n_eqn; i++){
                   yout[i]  = 0;
                   ypout[i] = 0;
              }
              g[1]   = 1.0;
              rho[1] = 1.0;
              hi = tout - x;
              ki = kold + 1;
      
              // Initialize w[*] for computing g[*]
              for (int i=0; i<ki; i++){
                  temp1 = i;
                  w[i+1] = 1.0/temp1;
              }

              // Compute g[*]
              term = 0.0;
              for (int j=1; j<ki; j++){
                  psijm1 = psi_[j];
                  gamma = (hi + term)/psijm1;
                  eta = hi/psijm1;
                  for (int i=0; i<ki+1-j; i++){
                    w[i+1] = gamma*w[i+1] - eta*w[i+2];
                  }
                  g[j+1] = w[1];
                  rho[j+1] = gamma*rho[j];
                  term = psijm1;
              }
      
              // Interpolate for the solution yout and for
              // the derivative of the solution ypout      
              for (int j=0; j<ki; j++){
                  i = ki+1-j;
                  for (int l = 0; l < n_eqn; l++) {
                      yout[l]  = yout[l]  + g[i+1]   * phi[l][i+1];
                      ypout[l] = ypout[l] + rho[i+1] * phi[l][i+1];
                  }
              }  
            
              for(int num=0; num<n_eqn; num++){
                  yout[num] = y[num] + hi*yout[num];
              }

              for(int num1=0; num1<n_eqn; num1++){
                  y[num1] = yout[num1];
              }
              
              State_    = DE_STATE.DE_DONE; // Set return code
              t         = tout;             // Set independent variable
              told      = t;                // Store independent variable
              OldPermit = PermitTOUT;
              return;                       // Normal exit
          }                        
  
          // If cannot go past output point and sufficiently close,
          // extrapolate and return
          if ( !PermitTOUT && ( fabs(tout-x) < fouru*fabs(x) ) ){
              h = tout - x;
              Accel(x, yy, yp);              // Compute derivative yp(x) -------------------------------------------------------
              for(int n1=0; n1<n_eqn; n1++){ // Extrapolate vector from x to tout
                  y[n1] = yy[n1] + h*yp[n1];
              }
              State_    = DE_STATE.DE_DONE; // Set return code
              t         = tout;             // Set independent variable
              told      = t;                // Store independent variable
              OldPermit = PermitTOUT;
              return;                       // Normal exit
          }
  
        // Test for too much work
        //   if (nostep >= maxnum)
        //       State_ = DE_STATE.DE_NUMSTEPS; % Too many steps
        //       if (stiff) 
        //           State_ = DE_STATE.DE_STIFF;% Stiffness suspected
        //       end
        //       y         = yy;                % Copy last step
        //       t         = x;
        //       told      = t;
        //       OldPermit = true;
        //       return;                        % Weak failure exit
        //   end
  
          // Limit step size, set weight vector and take a step
          
          if(fabs(h) > fabs(tend-x)){
              aux1 = fabs(tend-x);
          }else{
              aux1 = fabs(h);
          }
          h  = sign_(aux1, h);
          for (int l=0; l<n_eqn; l++){
              wt[l] = releps*fabs(yy[l]) + abseps;
          }
  
        //   Step
        //                                                                   
        // Begin block 0                                                     
        //                                                                   
        // Check if step size or error tolerance is too small for machine    
        // precision.  If first step, initialize phi array and estimate a    
        // starting step size. If step size is too small, determine an       
        // acceptable one.                                                   
        //                                                                   

        if (fabs(h) < fouru*fabs(x)){
            h = sign_(fouru*fabs(x),h);
            crash = true;
            return;           // Exit 
        }

        p5eps  = 0.5*epsilon;
        crash  = false;
        g[1]   = 1.0;
        g[2]   = 0.5;
        sig[1] = 1.0;

        ifail = 0;

        // If error tolerance is too small, increase it to an 
        // acceptable value.                                  

        round = 0.0;
        for (int l=0; l<n_eqn; l++){
            round = round + (y[l]*y[l])/(wt[l]*wt[l]);
        }

        round = twou*sqrt(round);
        if (p5eps<round){
            epsilon = 2.0*round*(1.0+fouru);
            crash = true;
            return;
        }

        if (start){
           // Initialize. Compute appropriate step size for first step. 
           Accel(x, y, yp); //--------------------------------------------------------------------------------------------
           sum = 0.0;
           for (int l=0; l<n_eqn; l++){
              phi[l][2] = yp[l];
              phi[l][3] = 0.0;
              sum = sum + (yp[l]*yp[l])/(wt[l]*wt[l]);
           }
             
           sum  = sqrt(sum);
           absh = fabs(h);
           if (epsilon < 16.0*sum*h*h){
              absh = 0.25*sqrt(epsilon/sum);
           }
           
           if(absh > fouru*fabs(x)){
                aux2 = absh;
           }else{
                aux2 = fouru*fabs(x);
           }

           h    = sign_(aux2, h);
           hold = 0.0;
           hnew = 0.0;
           k    = 1;
           kold = 0;
           start  = false;
           phase1 = true;
           nornd  = true;
           if (p5eps<=100.0*round){
              nornd = false;
              for (int l=0; l<n_eqn; l++){
                  phi[l][15]=0.0;
              }
           }
           return;
        }
                                                                   
        // End block 0 

        // Repeat blocks 1, 2 (and 3) until step is successful                                                                                
        while(true)
                                                                  
          // Begin block 1                                                   
          //                                                                 
          // Compute coefficients of formulas for this step. Avoid computing 
          // those quantities not changed when step size is not changed.                                                                      
  
          kp1 = k+1;
          kp2 = k+2;
          km1 = k-1;
          km2 = k-2;
  
          // ns is the number of steps taken with size h, including the 
          // current one. When k<ns, no coefficients change.           
  
          if (h != hold){
              ns = 0;
          }
          if (ns <= kold){
              ns = ns+1;
          }
          nsp1 = ns+1;
  
          if (k >= ns){
              // Compute those components of alpha[*],beta[*],psi[*],sig[*] 
              // which are changed                                          
              beta[ns] = 1.0;
              realns = ns;
              alpha[ns] = 1.0/realns;
              temp1 = h*realns;
              sig[nsp1] = 1.0;
              if (k >= nsp1)
                  for (int i=nsp1; i<k; i++){
                      im1   = i-1;
                      temp2 = psi_[im1];
                      psi_[im1] = temp1;
                      beta[i]  = beta[im1]*psi_[im1]/temp2;
                      temp1    = temp2 + h;
                      alpha[i] = h/temp1;
                      reali = i;
                      sig[i+1] = reali*alpha[i]*sig[i];
                  }
                      
              psi_[k] = temp1;
      
              // Compute coefficients g[*]; initialize v[*] and set w[*].
              if (ns>1){
                 // If order was raised, update diagonal part of v[*]
                  if (k>kold){
                      temp4 = k*kp1;
                      v[k] = 1.0/temp4;
                      nsm2 = ns-2;
                      for (int j=1; j<nsm2; j++){ 
                           i = k-j;
                           v[i] = v[i] - alpha[j+1]*v[i+1];
                      }
                  }
                      
                  // Update V[*] and set W[*]
                  limit1 = kp1 - ns;
                  temp5  = alpha[ns];
                  for (int iq=0; iq<limit1; iq++){
                      v[iq+1] = v[iq+1] - temp5*v[iq+2];
                      w[iq+1] = v[iq+1];
                  }
                  g[nsp1] = w[1];

              }else{
                 for (int iq=0; iq<k; iq++){
                     temp3 = iq*(iq+1);
                     v[iq+1] = 1.0/temp3;
                     w[iq+1] = v[iq+1];
                 }
              }
      
              // Compute the g[*] in the work vector w[*]
              nsp2 = ns + 2;
              if (kp1 >= nsp2){
                for (int i=nsp2; i<kp1; i++){
                    limit2 = kp2 - i;
                    temp6  = alpha[i-1];
                    for (int iq=0; iq<limit2; iq++){
                         w[iq+1] = w[iq+1] - temp6*w[iq+2];
                    }
                    g[i] = w[1];
                }
              }
          }

          // if K>=NS
  
          // End block 1
  
          // Begin block 2
  
          // Predict a solution p[*], evaluate derivatives using predicted
          // solution, estimate local error at order k and errors at orders
          // k, k-1, k-2 as if constant step size were used.  
  
          // Change phi to phi star
          if (k >= nsp1){
             for (int i=nsp1; i<k; i++){
                 temp1 = beta[i];
                 for (int l=0; l<n_eqn; l++){
                     phi[l][i] = temp1 * phi[l][i];
                 }
             }
          }
  
          // Predict solution and differences 
          for (int l=0; l<n_eqn; l++){
              phi[l][kp2] = phi[l][kp1];
              phi[l][kp1] = 0.0;
              p[l]        = 0.0;
          }
          for (int j=1; j<k; j++){
              i     = kp1 - j;
              ip1   = i+1;
              temp2 = g[i];
              for (int l=0; l<n_eqn; l++){
                  p[l]     = p[l] + temp2*phi[l][i];
                  phi[l][i] = phi[l][i] + phi[l][ip1];
              }
          }
          if (nornd){
              for(int jj=0; jj<n_eqn; jj++){
                    p[jj] = y[jj] + h*p[jj];
              }
          }else{
              for (int l=0; l<n_eqn; l++){
                  tau = h*p[l] - phi[l][15];
                  p[l] = y[l] + tau;
                  phi[l][16] = (p[l] - y[l]) - tau;
              }
          }

          xold = x;
          x = x + h;
          absh = fabs(h);
          Accel(x, p, yp); //------------------------------------------------------------------------------------------
  
          // Estimate errors at orders k, k-1, k-2 
          erkm2 = 0.0;
          erkm1 = 0.0;
          erk = 0.0;
  
          for (int l=0; l<n_eqn; l++){
              temp3 = 1.0/wt[l];
              temp4 = yp[l] - phi[l][1];
              if (km2 > 0){
                 erkm2 = erkm2 + ((phi[l][km1]+temp4)*temp3)*((phi[l][km1]+temp4)*temp3);
              }
              if (km2 >= 0){
                 erkm1 = erkm1 + ((phi[l][k]+temp4)*temp3)*((phi[l][k]+temp4)*temp3);
              }
              erk = erk + (temp4*temp3)*(temp4*temp3);
          }
  
          if (km2> 0){
            erkm2 = absh*sig[km1]*gstr[km2]*sqrt(erkm2);
          }
          if (km2>=0){
            erkm1 = absh*sig[k]*gstr[km1]*sqrt(erkm1);
          }
  
          temp5 = absh*sqrt(erk);
          err = temp5*(g[k]-g[kp1]);
          erk = temp5*sig[kp1]*gstr[k];
          knew = k;
  
          // Test if order should be lowered 
          if (km2 > 0){
              if(erkm1 > erkm2){
                    aux3 = erkm1;
              }else{
                    aux3 = erkm2;
              }
              if (aux3 <= erk){
                  knew = km1;
              }
          }
          if (km2 == 0){
              if (erkm1<=0.5*erk){
                  knew = km1;
              }
          }
  
          // End block 2
  
          // If step is successful continue with block 4, otherwise repeat
          // blocks 1 and 2 after executing block 3
          if(err <= epsilon){
                success = true;
          }else{
                success = false;
          }
  
          if (!success){
                // Begin block 3
    
                // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
                // 3rd consecutive failure, set order to 1. If step fails more
                // than 3 times, consider an optimal step size. Double error
                // tolerance and return if estimated step size is too small
                // for machine precision.
    
                // Restore x, phi[*,*] and psi[*]
                phase1 = false; 
                x = xold;
                for (int i=0; i<k; i++){
                    temp1 = 1.0/beta[i+1];
                    ip1 = i+1;
                    for (int l=0; l<n_eqn; l++){
                         phi[l][i+1]=temp1*(phi[l][i+1]-phi[l][ip1+1]);
                    }
                }
    
                if (k >= 2){
                    for (int i=1; i<k; i++){
                         psi_[i] = psi_[i+1] - h;
                    }
                }
    
                // On third failure, set order to one. 
                // Thereafter, use optimal step size   
                ifail = ifail + 1;
                temp2 = 0.5;
                if (ifail > 3){
                    if (p5eps < 0.25*erk){
                        temp2 = sqrt(p5eps/erk);
                    }
                }
                if (ifail >= 3){
                    knew = 1;
                }
                h = temp2*h;
                k = knew;
                if (fabs(h) < fouru*fabs(x)){
                    crash = true;
                    h = sign_(fouru*fabs(x), h);
                    epsilon = epsilon*2.0;
                    return;         // Exit 
                }
    
                // End block 3, return to start of block 1
    
                // end if(success)
  
              if (success){
                  break;
              }
          }
  


        // Begin block 4

        // The step is successful. Correct the predicted solution, evaluate
        // the derivatives using the corrected solution and update the
        // differences. Determine best order and step size for next step.

        kold = k;
        hold = h;

        // Correct and evaluate
        temp1 = h*g[kp1];
        if (nornd){
            for (int l=0; l<n_eqn; l++){
                y[l] = p[l] + temp1*(yp[l] - phi[l][1]);
            }
        }else{
            for (int l=0; l<n_eqn; l++){
                rho[l] = temp1*(yp[l] - phi[l][1]) - phi[l][16];
                y[l] = p[l] + rho[l];
                phi[l][15] = (y[l] - p[l]) - rho[l];
            }
        }
        Accel(x, y, yp); //-----------------------------------------------------------------------------------------------------

        // Update differences for next step 
        for (int l=0; l<n_eqn; l++){
            phi[l][kp1] = yp[l] - phi[l][2];
            phi[l][kp2] = phi[l][kp1] - phi[l][kp2];
        }
        for (int i=0; i<k; i++){
            for (int l=0; l<n_eqn; l++){
                phi[l][i+1] = phi[l][i+1] + phi[l][kp1];
            }
        }

        // Estimate error at order k+1 unless               
        // - in first phase when always raise order,        
        // - already decided to lower order,                
        // - step size not constant so estimate unreliable  
        erkp1 = 0.0;
        if ( (knew == km1) || (k == 12) ){
            phase1 = false;
        }

        if (phase1){ 
            k = kp1;
            erk = erkp1;
        }else{
            if (knew == km1){
                // lower order 
                k = km1;
                erk = erkm1;
            }else{
                if (kp1<=ns){
                    for (int l=0; l<n_eqn; l++){
                        erkp1 = erkp1 + (phi[l][kp2]/wt[l])*(phi[l][kp2]/wt[l]);
                    }
                    erkp1 = absh*gstr[kp1]*sqrt(erkp1);
                    // Using estimated error at order k+1, determine 
                    // appropriate order for next step               
                    if (k > 1){
                        if(erk > erkp1){
                            aux4 = erkp1;
                        }else{
                            aux4 = erk;
                        }
                        if (erkm1 <= aux4){
                            // lower order
                            k = km1; 
                            erk = erkm1;
                        }else if ((erkp1 < erk) && (k != 12)){
                            // raise order 
                            k = kp1;
                            erk = erkp1;
                        }
                    }else if (erkp1 < 0.5*erk){
                        // raise order
                        // Here erkp1 < erk < max(erkm1,ermk2) else    
                        // order would have been lowered in block 2.   
                        // Thus order is to be raised                  
                        k = kp1;
                        erk = erkp1;
                    } 
                }    
            }     
        }
            

        // With new order determine appropriate step size for next step
        if ( phase1 || (p5eps >= erk*two[k+1]) ){
             hnew = 2.0*h;
        }else{
            if (p5eps < erk){
                temp2 = k+1;
                r = p5eps/pow(erk,(1.0/temp2));
                if(0.9 > r){
                    aux5 = r; 
                }else{
                    aux5 = 0.9;
                }
                if(0.5 > aux5){
                    aux6 = 0.5; 
                }else{
                    aux6 = aux5;
                }
                hnew = absh*aux6;
                if(hnew > fouru*fabs(x)){
                    aux7 = hnew;
                }else{
                    aux7 = fouru*fabs(x);
                }
                hnew = sign_(aux7, h);
            }else{
                hnew = h;
            }
        }
        h = hnew; 

        // End block 4

          // Test for too small tolerances
          if (crash){
              State_    = DE_STATE.DE_BADACC;
              relerr    = epsilon*releps;       // Modify relative and absolute
              abserr    = epsilon*abseps;       // accuracy requirements
              for(int i=0; i<n_eqn; i++){
                  y[i] = yy[i];       
              }
              t         = x;
              told      = t;
              OldPermit = true;
              return; 
          }
                                    // Weak failure exit
   } // End step loop
  
    //   if ( State_==DE_STATE.DE_INVPARAM )
    //       error ('invalid parameters in DEInteg');
    //       exit; 
    //   end
    //   if ( State_==DE_STATE.DE_BADACC )
    //       warning ('on','Accuracy requirement not achieved in DEInteg');
    //   end
    //   if ( State_==DE_STATE.DE_STIFF )
    //       warning ('on','Stiff problem suspected in DEInteg');
    //   end
    //   if ( State_ >= DE_STATE.DE_DONE )
    //       break;
    //   end
    //   
    // end
}