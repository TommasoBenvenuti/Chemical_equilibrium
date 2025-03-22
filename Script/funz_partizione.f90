program funz_partizione
  USE momento_inerzia
  USE part_vib
  implicit none
  integer                       :: n_atomi = 14, n_vib
  integer                       :: i, j
  real                          :: k_b, T, sigma_a,sigma_b, pi, h, Eel_A, Eel_B, c, EvibA, EvibB, delta_E0, K_p
  real                          :: X_cmassa, Y_cmassa, Z_cmassa, M, In(3,3), det, lnf_rot, f_vib, funz_totA, funz_totB
  real, allocatable             :: X(:), Y(:), Z(:), freq_vib(:), Theta_vib(:)
  character, allocatable :: atomi(:)

  !==== Costanti in doppia precisione ====
  h    = 6.626068E-34! Costante di Planck in J·s
  k_b  = 1.380649E-23 ! Costante di Boltzmann in J/K
  T    = 298.15       ! Temperatura in Kelvin
  c    = 2.997925E10  ! Velocità della luce in cm/s

  ! ==== Allocazione array ====
  allocate(atomi(n_atomi), X(n_atomi), Y(n_atomi), Z(n_atomi))
  n_vib = 3 * n_atomi - 6
  allocate(freq_vib(n_vib), Theta_vib(n_vib))

  ! ==== Lettura primo file di dati (ag1.dat) ====
  open(unit=10, file='ag1.dat', status='old', action='read', iostat=i)
  if (i /= 0) then
    print *, 'Errore: impossibile aprire il file ag1.dat'
    stop
  end if

  read(10,*)
  do i = 1, n_atomi
     read(10,*) atomi(i), X(i), Y(i), Z(i)
  end do
  read(10,*) ! Ignora una riga
  read(10,*) ! Ignora una seconda riga
  do j = 1, n_vib
     read(10,*) freq_vib(j)
  end do
  read(10,*) Eel_A
  read(10,*) sigma_a
  close(10)

  Eel_A = Eel_A * 4.3597E-18

  ! ==== Calcolo e chiamata alle funzioni ====
  call MASSA(atomi, M)
  call coord_CM(atomi, M, X, X_cmassa, Y, Y_cmassa, Z, Z_cmassa)
  call calcola_tensore(atomi, M, X, Y, Z, X_cmassa, Y_cmassa, Z_cmassa, In)
  call calcola_det(In, det)
  call calcola_partizione_vib(n_vib, freq_vib, Theta_vib, f_vib, EvibA)

  funz_totA = f_vib * sqrt(det/sigma_a)
  write(*,*) 'funz_totA =', funz_totA

  ! ==== Lettura secondo file di dati (ag2.dat) ====
  open(unit=10, file='ag2.dat', status='old', action='read', iostat=i)
  if (i /= 0) then
    print *, 'Errore: impossibile aprire il file ag2.dat'
    stop
  end if

  read(10,*)
  DO i = 1, n_atomi
     read(10,*) atomi(i), X(i), Y(i), Z(i)
  ENDDO
  read(10,*)
  read(10,*)
  DO j = 1, n_vib
     read(10,*) freq_vib(j)
  ENDDO
  read(10,*) Eel_B
  read(10,*) sigma_b
  close(10)

  Eel_B = Eel_B * 4.3597E-18

  ! ==== Calcolo e chiamata funzioni per seconda molecola ====
  call MASSA(atomi, M)
  call coord_CM(atomi, M, X, X_cmassa, Y, Y_cmassa, Z, Z_cmassa)
  call calcola_tensore(atomi, M, X, Y, Z, X_cmassa, Y_cmassa, Z_cmassa, In)
  call calcola_det(In, det)
  call calcola_partizione_vib(n_vib, freq_vib, Theta_vib, f_vib, EvibB)

  funz_totB = f_vib * sqrt(det/sigma_b) 
  write(*,*) 'funz_totB =', funz_totB

  ! ==== Calcolo della costante di equilibrio ====
  delta_E0 = (EvibB + Eel_B - EvibA - Eel_A) / (k_b * T)
  K_p = exp(-delta_E0) * funz_totB / funz_totA
  write(*,*) 'kp =', K_p

  ! ==== Scrittura risultati in un file ====
  open(unit=20, file='risultati.txt', status='replace', action='write')
  write(20,*) '==========================================================='
  write(20,*) '            RISULTATI DEL CALCOLO ENERGETICO              '
  write(20,*) '==========================================================='
  write(20,*) '--------------------- Molecola 1 --------------------------'
  write(20, '(A, E12.3, A,A, E12.3, A,A, E12.3)') 'EvibA (J): ', EvibA, '---', &
 'funz_totA: ', funz_totA, '---', 'Eel_A (J): ', Eel_A

  write(20,*) '--------------------- Molecola 2 --------------------------'
  write(20, '(A, E12.3,A, A, E12.3, A, A, E12.3)') 'EvibB (J): ', EvibB, '---', &
  'funz_totB: ', funz_totB, '---', 'Eel_B (J): ', Eel_B

  write(20,*) '==========================================================='
  write(20,*) '                  Costante Equilibrio                     '
  write(20,*) '==========================================================='
  write(20,'(A, E12.5)') 'Costante di Equilibrio: ', K_p
  write(20,*) '==========================================================='
  close(20)

  ! Deallocazione memoria
  deallocate(atomi, X, Y, Z, freq_vib, Theta_vib)

end program funz_partizione
