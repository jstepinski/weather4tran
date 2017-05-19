MODULE funcs

! Template for the calculator functions used in the Matsuno scheme

INTERFACE

            FUNCTION pi_a_calc(pi_a_prev)
              REAL, INTENT(IN) :: pi_a_prev(:,:)
              REAL pi_a_calc(size(pi_a_prev,1), size(pi_a_prev,2))
            END FUNCTION
            
            FUNCTION sigmaDot_calc(pi_a_prev, pi_a_curr, n)
              REAL, INTENT(IN) :: pi_a_prev(:,:)
              REAL, INTENT(IN) :: pi_a_curr(:,:)
              REAL sigmaDot_calc(size(pi_a_prev,1), size(pi_a_prev,2), n)
            END FUNCTION
            
            FUNCTION theta_interp_calc(theta)
              REAL, INTENT(IN) :: theta(:,:,:)
              REAL theta_interp_calc(size(theta,1), size(theta,2), size(theta,3)+1)
            END FUNCTION
            
            FUNCTION theta_calc(pi_a_prev, pi_a_curr, theta_prev, theta_curr, theta_interp, sigmaDot)
              REAL, INTENT(IN) :: pi_a_prev(:,:)
              REAL, INTENT(IN) :: pi_a_curr(:,:)
              REAL, INTENT(IN) :: theta_prev(:,:,:)
              REAL, INTENT(IN) :: theta_curr(:,:,:)
              REAL, INTENT(IN) :: theta_interp(:,:,:)
              REAL, INTENT(IN) :: sigmaDot(:,:,:)
              
              REAL theta_calc(size(theta_prev,1), size(theta_prev,2), size(theta_prev,3))
            END FUNCTION
            
            FUNCTION phi_calc(theta)
              REAL, INTENT(IN) :: theta(:,:,:)
              REAL phi_calc(size(theta,1), size(theta,2), size(theta,3))
            END FUNCTION
            
            FUNCTION uwind_calc(uwind_prev, uwind_curr, vwind, pi_a_prev, pi_a_prev2, pi_a_curr, pi_a_est, sigmaDot, phi, phi_prev, theta)
              REAL, INTENT(IN) :: uwind_prev(:,:,:)
              REAL, INTENT(IN) :: uwind_curr(:,:,:)
              REAL, INTENT(IN) :: vwind(:,:,:)
              REAL, INTENT(IN) :: pi_a_prev(:,:)
              REAL, INTENT(IN) :: pi_a_prev2(:,:)
              REAL, INTENT(IN) :: pi_a_curr(:,:)
              REAL, INTENT(IN) :: pi_a_est(:,:)
              REAL, INTENT(IN) :: sigmaDot(:,:,:)
              REAL, INTENT(IN) :: phi(:,:,:)
              REAL, INTENT(IN) :: phi_prev(:,:,:)
              REAL, INTENT(IN) :: theta(:,:,:)
              
              REAL uwind_calc(size(uwind_prev,1), size(uwind_prev,2), size(uwind_prev,3))
            END FUNCTION
            
            FUNCTION vwind_calc(vwind_prev, vwind_curr, uwind, pi_a_prev, pi_a_prev2, pi_a_curr, pi_a_est, sigmaDot, phi, phi_prev, theta)
              REAL, INTENT(IN) :: vwind_prev(:,:,:)
              REAL, INTENT(IN) :: vwind_curr(:,:,:)
              REAL, INTENT(IN) :: uwind(:,:,:)
              REAL,  INTENT(IN) :: pi_a_prev(:,:)
              REAL,  INTENT(IN) :: pi_a_prev2(:,:)
              REAL,  INTENT(IN) :: pi_a_curr(:,:)
              REAL,  INTENT(IN) :: pi_a_est(:,:)
              REAL,  INTENT(IN) :: sigmaDot(:,:,:)
              REAL,  INTENT(IN) :: phi(:,:,:)
              REAL,  INTENT(IN) :: phi_prev(:,:,:)
              REAL,  INTENT(IN) :: theta(:,:,:)
              
              REAL vwind_calc(size(vwind_prev,1), size(vwind_prev,2), size(vwind_prev,3))
            END FUNCTION

        END INTERFACE

END MODULE
