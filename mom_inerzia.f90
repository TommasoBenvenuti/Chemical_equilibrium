MODULE momento_inerzia
  implicit none
CONTAINS
  
  SUBROUTINE MASSA(atomi,M)
    IMPLICIT NONE
    character, INTENT(IN), DIMENSION(:) :: atomi
    real, INTENT(OUT) :: M 
    integer :: n_atomi=14,mX
    integer :: i

   ! ==== Calcolo la massa totale della molecola ====

    M=0
    do i=1,n_atomi
       if (atomi(i)=='C') then
          mX=12
       elseif (atomi(i)=='H') then
          mX=1
       else
          mX=16
       endif
       M=M+mX
    enddo
  END SUBROUTINE MASSA
  
  SUBROUTINE coord_CM(atomi,M,X,X_cmassa,Y,Y_cmassa,Z,Z_cmassa)
    implicit none
    character, INTENT(IN), DIMENSION(:) :: atomi(:)
    real, INTENT(IN), DIMENSION(:) :: X(:),Y(:),Z(:)
    real, INTENT(IN) :: M
    real, INTENT(OUT) :: X_cmassa,Y_cmassa,Z_cmassa
    real :: mX
    integer :: n_atomi=14
    integer :: i

   ! ===== Calcolo la posizione del CM ====

    X_cmassa=0
    Y_cmassa=0
    Z_cmassa=0

    DO i=1,n_atomi
       if (atomi(i)=='C') then
          mX=12
       elseif (atomi(i)=='H') then
          mX=1
       else
          mX=16
       endif
       X_cmassa=X_cmassa+mX*X(i)
       Y_cmassa=Y_cmassa+mX*Y(i)
       Z_cmassa=Z_cmassa+mX*Z(i)
    ENDDO

    X_cmassa=X_cmassa/M
    Y_cmassa=Y_cmassa/M
    Z_cmassa=Z_cmassa/M   
  END SUBROUTINE COORD_CM
  
  SUBROUTINE calcola_tensore(atomi,M,X,Y,Z,X_cmassa,Y_cmassa,Z_cmassa,In)
    implicit none
    integer :: i,j,s,n_atomi=14
    character, INTENT(IN), DIMENSION(:) :: atomi(:)
    real, Intent(IN), DIMENSION(:) :: X(:),Y(:),Z(:)
    real, INTENT(IN) :: M,X_cmassa,Y_cmassa,Z_cmassa
    real, INTENT(OUT) :: In(3,3)
    real, DIMENSION(14) :: xS(14),yS(14),zS(14)
    real :: A(3,14)
    real :: mX

! ==== Calcolo le posizioni traslate degli atomi ====

    DO i=1,n_atomi
       xS(i)=(X(i)-X_cmassa)
       yS(i)=(Y(i)-Y_cmassa)
       zS(i)=(Z(i)-Z_cmassa)
    ENDDO
! In una matrice 3*14 stampo le coordinate traslate degli atomi    
    DO i=1,3
       DO j=1,14
          A(i,j)=0.
       ENDDO
    ENDDO
   DO j=1,14
             A(1,j)=xS(j)
             A(2,j)=yS(j)
             A(3,j)=zS(j)
   ENDDO
! ==== Calcolo tensore di inerzia ====    
    DO i=1,3
       DO j=1,3
          In(i,j)=0.
       ENDDO
    ENDDO
    DO i=1,3
       DO j=1,3
         DO s=1,n_atomi
             if (atomi(s)=='C') then
                mX=12
             elseif (atomi(s)=='H') then
                mX=1
             else
                mX=16
             endif
             if (i/=j) then
                In(i,j)=In(i,j)-A(i,s)*A(j,s)*mX
             else
                In(1,1) = In(1,1) + mX * (yS(s)**2 + zS(s)**2)
                In(2,2) = In(2,2) + mX * (xS(s)**2 + zS(s)**2)
                In(3,3) = In(3,3) + mX * (xS(s)**2 + yS(s)**2)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    write(*,*) 'In(i,j)'
    do i=1,3
       write(*, '(3F10.3)') (In(i, j), j = 1, 3)
    enddo
  END SUBROUTINE calcola_tensore

  SUBROUTINE calcola_det(In,det)
    real, INTENT(IN) :: In(3,3)
    real, INTENT(OUT) :: det
    det=0.
    det = In(1,1)*In(2,2)*In(3,3) - In(1,2)*In(2,3)*In(3,1) - In(1,3)*In(2,1)*In(3,2)  &
      - In(1,1)*In(2,3)*In(3,2) - In(1,2)*In(2,1)*In(3,3) - In(1,3)*In(2,2)*In(3,1)
    write(*,*) 'det=',det
  END SUBROUTINE calcola_det
END MODULE momento_inerzia

