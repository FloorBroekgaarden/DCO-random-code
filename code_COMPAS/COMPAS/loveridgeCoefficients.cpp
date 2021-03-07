//
//  loveridgeCoefficients.cpp


#include "loveridgeCoefficients.h"

/*
 !*********************************************************************************************************************************
 !> \brief This function computes log[BE/erg] as a function of log[Z], Mzams, M, log[R/Ro] and GB.
 !!
 !! \param logZ   10-base log of metallicity (solar = log[0.02]);
 !! \param Mzams  stellar ZAMS mass (in solar masses);
 !! \param M      current stellar mass (in solar masses);
 !! \param logR   10-base log of stellar radius, expressed in solar radii (log[R/Ro]);
 !! \param GB     giant branch: 1: RGB, 2: AGB (CO core exists).
 !! \param ignore_massloss  if true, ignore the factor Lambda that corrects for mass loss (logical, optional; default: false)
 !! \retval calc_logBE  The 10-base log of the absolute value of the envelope binding energy expressed in erg.
 
 function calc_logBE(logZ, Mzams, M, logR, GB,  ignore_massloss)
 use kinds
 use BE_data
 
 implicit none
 integer, intent(in) :: GB
 real(double), intent(in) :: logZ,Mzams,M,logR
 logical, optional, intent(in) :: ignore_massloss
 integer :: ii,iz,ig
 real(double) :: calc_logBE,logM,dz,dzi,logRbound,dlogBE,lambda
 
 logM = log10(M)
 iz = 0
 
 ! Find the nearest logZ in the grid:
 dz = huge(dz)
 do ii=1,nz
 dzi = abs(log10(zs(ii))-logZ)
 if(dzi.lt.dz) then
 dz = dzi
 iz = ii
 end if
 end do
 
 
 
 ! Find the group from the mass and evolutionary state (ig = 1-4):
 if(logM.le.LMHMbs(iz)) then                           ! group LM (low mass)
 if(GB.eq.1) then                                   ! group LMR (low-mass RGB)
 ig = 1                                          ! group LMR1 (low-mass RGB 1)
 
 ! Compute boundary radius between LMR1 an LMR2:
 logRbound = 0.0_dbl
 do ii=0,nRGBbc-1
 logRbound = logRbound + RGBb_coef(iz,ii)*logM**ii
 end do
 if(logR.gt.logRbound) ig = 2                    ! group LMR2 (low-mass RGB 2)
 else                                               ! group LMA (low-mass AGB)
 ig = 3
 end if
 else                                                  ! group HM (high mass)
 ig = 4
 end if
 
 
 ! Compute the 10-logarithm of the binding energy:
 calc_logBE = 0.0_dbl
 do ii=1,ndat(iz,ig)
 dlogBE = alphas(iz,ig,ii) * (logM+tiny(logM))**dble(ms(iz,ig,ii)) * (logR+tiny(logR))**dble(rs(iz,ig,ii))  ! Avoid 0**0
 calc_logBE = calc_logBE + dlogBE
 end do
 
 
 ! Compute and apply the mass-loss correction factor Lambda:
 lambda = 1.0_dbl + 0.25_dbl * ((Mzams-M)/Mzams)**2
 if(present(ignore_massloss)) then
 if(ignore_massloss) lambda = 1.0_dbl  ! If the optional parameter ignore_massloss is present and true, set Lambda=1
 end if
 
 calc_logBE = lambda * calc_logBE
 
 
 ! BE was originally expressed in erg/solar mass to avoid large numbers, so we need to convert to erg here:
 calc_logBE = calc_logBE + logBE0
 
 end function calc_logBE
 !*********************************************************************************************************************************
 


*/