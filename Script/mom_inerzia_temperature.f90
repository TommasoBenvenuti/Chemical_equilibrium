MODULE momento_inerzia
  implicit none
CONTAINS

  SUBROUTINE ASSEGNAMASSA(simbolo, massa)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: simbolo
    REAL, INTENT(OUT) :: massa
    
    SELECT CASE(TRIM(simbolo))
      CASE ('H')
        massa = 1.008
      CASE ('C')
        massa = 12.011
      CASE ('N')
        massa = 14.007
      CASE ('O')
        massa = 16.000
      CASE ('F')
        massa = 18.998
      CASE ('P')
        massa = 30.974
      CASE ('S')
        massa = 32.065
      CASE ('Cl')
        massa = 35.453
      CASE ('Fe')
        massa = 55.845
      CASE ('Co')
        massa = 58.933
      CASE ('Br')
        massa = 79.904
      CASE ('I')
        massa = 126.90
      CASE DEFAULT
        PRINT *, 'Attenzione: Simbolo atomico non riconosciuto -> ', simbolo
        massa = 0.0
    END SELECT
  END SUBROUTINE ASSEGNAMASSA

  SUBROUTINE MASSA(atomi, M)
    IMPLICIT NONE
    CHARACTER, INTENT(IN), DIMENSION(:) :: atomi
    REAL, INTENT(OUT) :: M 
    INTEGER :: i
    REAL :: massaX
    
    M = 0.0
    DO i = 1, SIZE(atomi)
       CALL ASSEGNAMASSA(atomi(i), massaX)
       M = M + massaX
    ENDDO
  END SUBROUTINE MASSA

  SUBROUTINE coord_CM(atomi, M, X, X_cmassa, Y, Y_cmassa, Z, Z_cmassa)
    IMPLICIT NONE
    CHARACTER, INTENT(IN), DIMENSION(:) :: atomi
    REAL, INTENT(IN), DIMENSION(:) :: X, Y, Z
    REAL, INTENT(IN) :: M
    REAL, INTENT(OUT) :: X_cmassa, Y_cmassa, Z_cmassa
    REAL :: massaX
    INTEGER :: i

    X_cmassa = 0.0
    Y_cmassa = 0.0
    Z_cmassa = 0.0

    DO i = 1, SIZE(atomi)
       CALL ASSEGNAMASSA(atomi(i), massaX)
       X_cmassa = X_cmassa + massaX * X(i)
       Y_cmassa = Y_cmassa + massaX * Y(i)
       Z_cmassa = Z_cmassa + massaX * Z(i)
    ENDDO

    X_cmassa = X_cmassa / M
    Y_cmassa = Y_cmassa / M
    Z_cmassa = Z_cmassa / M   
  END SUBROUTINE coord_CM

  SUBROUTINE calcola_tensore(atomi, M, X, Y, Z, X_cmassa, Y_cmassa, Z_cmassa, In)
    IMPLICIT NONE
    INTEGER :: i, j, s
    CHARACTER, INTENT(IN), DIMENSION(:) :: atomi
    REAL, INTENT(IN), DIMENSION(:) :: X, Y, Z
    REAL, INTENT(IN) :: M, X_cmassa, Y_cmassa, Z_cmassa
    REAL, INTENT(OUT) :: In(3,3)
    REAL, ALLOCATABLE :: xS(:), yS(:), zS(:)
    REAL :: massaX

    ALLOCATE(xS(SIZE(atomi)), yS(SIZE(atomi)), zS(SIZE(atomi)))

    DO i = 1, SIZE(atomi)
       xS(i) = X(i) - X_cmassa
       yS(i) = Y(i) - Y_cmassa
       zS(i) = Z(i) - Z_cmassa
    ENDDO

    In = 0.0
    DO s = 1, SIZE(atomi)
       CALL ASSEGNAMASSA(atomi(s), massaX)
       In(1,1) = In(1,1) + massaX * (yS(s)**2 + zS(s)**2)
       In(2,2) = In(2,2) + massaX * (xS(s)**2 + zS(s)**2)
       In(3,3) = In(3,3) + massaX * (xS(s)**2 + yS(s)**2)
       In(1,2) = In(1,2) - massaX * xS(s) * yS(s)
       In(1,3) = In(1,3) - massaX * xS(s) * zS(s)
       In(2,3) = In(2,3) - massaX * yS(s) * zS(s)
    ENDDO
    In(2,1) = In(1,2)
    In(3,1) = In(1,3)
    In(3,2) = In(2,3)

    DEALLOCATE(xS, yS, zS)
  END SUBROUTINE calcola_tensore

  SUBROUTINE calcola_det(In, det)
    REAL, INTENT(IN) :: In(3,3)
    REAL, INTENT(OUT) :: det

    det = In(1,1)*In(2,2)*In(3,3) - In(1,2)*In(2,3)*In(3,1) - &
          In(1,3)*In(2,1)*In(3,2) - In(1,1)*In(2,3)*In(3,2) - &
          In(1,2)*In(2,1)*In(3,3) - In(1,3)*In(2,2)*In(3,1)
  END SUBROUTINE calcola_det

END MODULE momento_inerzia
