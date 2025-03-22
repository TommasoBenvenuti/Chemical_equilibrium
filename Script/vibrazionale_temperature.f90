MODULE part_vib
  implicit none
CONTAINS
  SUBROUTINE calcola_partizione_vib(n_vib,freq_vib,Theta_vib,f_vib, Evib0, T)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n_vib, T
    REAL, INTENT(out)   :: Evib0
    REAL, INTENT(in), DIMENSION(:) :: freq_vib(:)
    REAL, INTENT(out), DIMENSION(:) :: Theta_vib
    REAL, INTENT(out), DIMENSION(:) :: f_vib(:)
    REAL :: h,pi,k_b,c
    INTEGER :: i, idT

  ! ==== Costanti ====
   
    h    = 6.626068E-34 ! Costante di Planck in J·s
    k_b  = 1.380649E-23 ! Costante di Boltzmann in J/K
    c    = 2.997925E10  ! Velocità della luce in cm/s

  ! ==== Calcolo delle temperature caratteristiche vibrazionali ====

    DO i=1, n_vib
       Theta_vib(i) = h*c*freq_vib(i)/(k_b)   
    ENDDO

  ! ==== Calcolo energia punto zero vibrazionale ====
  
 Evib0 = 0.0
    DO i=1,n_vib
     Evib0=Evib0+h*freq_vib(i)*c/(2)
 ENDDO
    
  ! ==== Calcolo della funzione di partizione vibrazionale ====

    f_vib = 1.0
    DO i=1, n_vib
       f_vib= f_vib * ( 1.0 / (1.0 - exp(-Theta_vib(i) / T )) )  
    ENDDO

    
    !write(*,*) Theta_vib
    !write(*,*) 'f_vib=',f_vib
    !write(*,*) freq_vib

  END SUBROUTINE
END MODULE
