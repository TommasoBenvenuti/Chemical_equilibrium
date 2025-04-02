program funz_partizione
  USE momento_inerzia
  USE part_vib
  implicit none
  integer                       :: n_atomi, n_vib
  integer                       :: i, j, idT, num_temperatures, T
  real                          :: k_b, sigma_a, sigma_b, pi, h, Eel_A, Eel_B, c, EvibA, EvibB, delta_E0, R
  real                          :: X_cmassa, Y_cmassa, Z_cmassa, M, In(3,3), det, lnf_rot
  real, allocatable             :: X(:), Y(:), Z(:), freq_vib(:), Theta_vib(:), funz_totA(:), funz_totB(:)
  real, allocatable             :: f_vib(:), K_p(:), deltaG(:)
  character, allocatable        :: atomi(:)
  character(len=100)            :: line
  
  !==== Costanti in doppia precisione ====
  h    = 6.626068E-34 ! Costante di Planck in J·s
  k_b  = 1.380649E-23 ! Costante di Boltzmann in J/K
  c    = 2.997925E10  ! Velocità della luce in cm/s
  R    = 8.31         ! J/(mol*K)
  num_temperatures = 80 
  
  ! ==== Determinazione del numero di atomi dal file ag1.dat ====
  open(unit=10, file='ag1.dat', status='old', action='read', iostat=i)
  read(10,*)
  if (i /= 0) then
    print *, 'Errore: impossibile aprire il file ag1.dat'
    stop
  end if

  n_atomi = 0
  do
    read(10, '(A)', iostat=i) line
    if (i /= 0) then
      print *, 'Errore: fine file prima di trovare "freq"'
      stop
    end if
    if (index(line, 'freq') > 0) exit
    n_atomi = n_atomi + 1
  end do
  rewind(10)

print *, n_atomi
  ! ==== Allocazione degli array ====
  allocate(atomi(n_atomi), X(n_atomi), Y(n_atomi), Z(n_atomi))
  n_vib = 3 * n_atomi - 6
  allocate(freq_vib(n_vib), Theta_vib(n_vib))

  ! ==== Lettura dei dati dal file ag1.dat ====
  read(10,*)
  do i = 1, n_atomi
    read(10,*) atomi(i), X(i), Y(i), Z(i)
  end do
  do
    read(10, '(A)') line
    if (index(line, 'freq') > 0) exit
  end do
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

  allocate(funz_totA(num_temperatures), f_vib(num_temperatures))
  DO idT = 1, num_temperatures
      T = 20 + (idT-1)*10
      call calcola_partizione_vib(n_vib, freq_vib, Theta_vib, f_vib, EvibA, T)
      funz_totA(idT) = f_vib(idT) * sqrt(det/sigma_a)
  ENDDO
  deallocate(f_vib)

  allocate(funz_totB(num_temperatures), f_vib(num_temperatures))

  ! ==== Lettura del file ag2.dat ====
  open(unit=10, file='ag2.dat', status='old', action='read', iostat=i)
  read(10,*)
  if (i /= 0) then
    print *, 'Errore: impossibile aprire il file ag2.dat'
    stop
  end if

  n_atomi = 0
  do
    read(10, '(A)', iostat=i) line
    if (i /= 0) then
      print *, 'Errore: fine file prima di trovare "freq"'
      stop
    end if
    if (index(line, 'freq') > 0) exit
    n_atomi = n_atomi + 1
  end do
  rewind(10)

  ! ==== Lettura dei dati dal file ag2.dat ====
  read(10,*)
  do i = 1, n_atomi
    read(10,*) atomi(i), X(i), Y(i), Z(i)
  end do
  do
    read(10, '(A)') line
    if (index(line, 'freq') > 0) exit
  end do
  do j = 1, n_vib
    read(10,*) freq_vib(j)
  end do
  read(10,*) Eel_B
  read(10,*) sigma_b
  close(10)

  Eel_B = Eel_B * 4.3597E-18

  ! ==== Calcolo per la seconda molecola ====
  call MASSA(atomi, M)
  call coord_CM(atomi, M, X, X_cmassa, Y, Y_cmassa, Z, Z_cmassa)
  call calcola_tensore(atomi, M, X, Y, Z, X_cmassa, Y_cmassa, Z_cmassa, In)
  call calcola_det(In, det)

  DO idT = 1, num_temperatures
      T = 20 + (idT-1)*10
      call calcola_partizione_vib(n_vib, freq_vib, Theta_vib, f_vib, EvibB, T)
      funz_totB(idT) = f_vib(idT) * sqrt(det/sigma_B)
  ENDDO

  ! ==== Apertura del file per scrivere i risultati ====
  open(unit=20, file='risultati_kp.dat', status='replace', action='write')
  allocate(K_p(num_temperatures), deltaG(num_temperatures))

  
  ! Intestazione
  write(20,*) '==========================================================='
  write(20,*) '            RISULTATI COSTANTE DI EQUILIBRIO               '
  write(20,*) '==========================================================='
  write(20,*) 'Temperatura (K)   Kp   deltaG'
  write(20,*) '-------------------------'

  DO idT = 1, num_temperatures
      T = 20 + (idT-1)*10
      K_p(idT) = exp(-(Eel_B + EvibB - EvibA - Eel_A) / (k_b * T)) * funz_totB(idT) / funz_totA(idT)
      deltaG(idT) = -R * T * log(K_p(idT))
      write(20, '(F10.3,2X,F10.3,2X,F10.3)') real(T), K_p(idT), deltaG(idT)
  ENDDO
  close(20)

  deallocate(atomi, X, Y, Z, freq_vib, Theta_vib, funz_totA, funz_totB, K_p, deltaG)
end program funz_partizione
