MODULE routines_mod
    IMPLICIT NONE

CONTAINS

    ! --------------------------------------------------------------------------
    ! DADOSMOD: Lê dados do arquivo DADOSMOD.txt
    ! --------------------------------------------------------------------------
    SUBROUTINE DADOSMOD(filename, n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI)
        CHARACTER(LEN=*), INTENT(IN)  :: filename
        INTEGER,          INTENT(OUT) :: n
        REAL(8), ALLOCATABLE, INTENT(OUT) :: COMPI(:), RAIOI(:), ESPHI(:), MEY(:)
        REAL(8),          INTENT(OUT) :: MI, RO, PI

        INTEGER :: iunit, i, ios
        REAL(8) :: t_MI, t_RO, t_PI

        iunit = 10
        OPEN(UNIT=iunit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            WRITE(*,*) 'ERRO: Nao foi possivel abrir o arquivo ', TRIM(filename)
            STOP
        END IF

        READ(iunit, *, IOSTAT=ios) n
        IF (ios /= 0) THEN
            WRITE(*,*) 'ERRO: Falha ao ler n na primeira linha de DADOSMOD.txt'
            STOP
        END IF

        ALLOCATE(COMPI(n), RAIOI(n), ESPHI(n), MEY(n))

        MI = 0.0D0
        RO = 0.0D0
        PI = 0.0D0

        DO i = 1, n
            READ(iunit, *, IOSTAT=ios) COMPI(i), RAIOI(i), ESPHI(i), MEY(i), &
                                       t_MI, t_RO, t_PI
            IF (ios /= 0) THEN
                WRITE(*,'(A,I4,A)') 'ERRO: Falha ao ler linha ', i+1, ' de DADOSMOD.txt'
                STOP
            END IF
            IF (i == 1) THEN
                MI = t_MI
                RO = t_RO
                PI = t_PI
            END IF
        END DO

        CLOSE(iunit)
    END SUBROUTINE DADOSMOD

    ! --------------------------------------------------------------------------
    ! GEOMAT: Calcula Resistência, Indutância e Compliância; salva GEOMAT.csv
    ! --------------------------------------------------------------------------
    SUBROUTINE GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI, &
                      RESISTI, INDUCTI, COMPLII)
        INTEGER, INTENT(IN)  :: n
        REAL(8), INTENT(IN)  :: COMPI(n), RAIOI(n), ESPHI(n), MEY(n)
        REAL(8), INTENT(IN)  :: MI, RO, PI
        REAL(8), INTENT(OUT) :: RESISTI(n), INDUCTI(n), COMPLII(n)

        INTEGER :: m, iunit
        REAL(8) :: raio, raio2, raio3, raio4

        DO m = 1, n
            raio  = RAIOI(m)
            raio2 = raio**2
            raio3 = raio**3
            raio4 = raio**4

            RESISTI(m) = (8.0D0 * MI * COMPI(m)) / (PI * raio4)
            INDUCTI(m) = (9.0D0 * RO * COMPI(m)) / (4.0D0 * PI * raio2)
            COMPLII(m) = (3.0D0 * PI * COMPI(m) * raio3) / (2.0D0 * MEY(m) * ESPHI(m))
        END DO

        ! Salvar CSV
        iunit = 20
        OPEN(UNIT=iunit, FILE='GEOMAT.csv', STATUS='REPLACE', ACTION='WRITE')
        WRITE(iunit,'(A)') 'Trecho,RESISTI,INDUCTI,COMPLII'
        DO m = 1, n
            WRITE(iunit,'(I6,3(",",ES15.6E2))') m, RESISTI(m), INDUCTI(m), COMPLII(m)
        END DO
        CLOSE(iunit)
    END SUBROUTINE GEOMAT

    ! --------------------------------------------------------------------------
    ! PULSOPQ: Calcula Pressão (PIN) e Vazão (QIN) de entrada no tempo IT*HT
    ! --------------------------------------------------------------------------
    SUBROUTINE PULSOPQ(IT, HT, PI, PIN, QIN)
        INTEGER, INTENT(IN)  :: IT
        REAL(8), INTENT(IN)  :: HT, PI
        REAL(8), INTENT(OUT) :: PIN, QIN

        REAL(8), PARAMETER :: TS   = 0.3D0
        REAL(8), PARAMETER :: TCC  = 1.0D0
        REAL(8), PARAMETER :: PS   = 120.0D0
        REAL(8), PARAMETER :: PD   = 80.0D0
        REAL(8), PARAMETER :: QMAX = 300.0D0

        REAL(8) :: TI, TI_MOD, A, B

        TI = IT * HT
        TI_MOD = MODULO(TI, TCC)  ! Torna o pulso periódico

        IF (TI_mod <= TS) THEN
            PIN = (PD + ((PS - PD) / TS) * TI_MOD) 
            QIN = QMAX * (DSIN(PI * TI_mod / TS))**2
        ELSE
            A = (PS - PD) / (TS - TCC)
            B = (PD * TS - PS * TCC) / (TS - TCC)
            PIN = (A * TI_mod + B) 
            QIN = 0.0D0
        END IF
    END SUBROUTINE PULSOPQ

    ! --------------------------------------------------------------------------
    ! FORMAB: Monta a matriz A e o vetor b do sistema
    ! --------------------------------------------------------------------------
    SUBROUTINE FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PIN, QIN, A, b)
        INTEGER, INTENT(IN)    :: n
        REAL(8), INTENT(IN)    :: RESISTI(n), INDUCTI(n), COMPLII(n)
        REAL(8), INTENT(IN)    :: HT, PIN, QIN
        REAL(8), INTENT(INOUT) :: A(2*n, 2*n), b(2*n)

        A(:,:) = 0.0D0
        b(:)   = 0.0D0

        
        ! --- LINHA 1 ---
        A(1,4)   =  1.0D0 / COMPLII(1)
        A(1,10)  =  1.0D0 / COMPLII(1)
        A(1,16)  =  1.0D0 / COMPLII(1)
        ! --- LINHA 2 ---
        A(2,1)   =  1.0D0 / INDUCTI(1)
        A(2,4)   = -RESISTI(1) / INDUCTI(1)
        A(2,10)  = -RESISTI(1) / INDUCTI(1)
        A(2,16)  = -RESISTI(1) / INDUCTI(1)
        ! --- LINHA 3 ---
        A(3,4)   =  1.0D0 / COMPLII(2)
        A(3,6)   = -1.0D0 / COMPLII(2)
        A(3,8)   = -1.0D0 / COMPLII(2)
        ! --- LINHA 4 ---
        A(4,1)   =  1.0D0 / INDUCTI(2)
        A(4,3)   = -1.0D0 / INDUCTI(2)
        A(4,4)   = -RESISTI(2) / INDUCTI(2)
        ! --- LINHA 5 ---
        A(5,6)   =  1.0D0 / COMPLII(3)
        ! --- LINHA 6 ---
        A(6,3)   =  1.0D0 / INDUCTI(3)
        A(6,5)   = -1.0D0 / INDUCTI(3)
        A(6,6)   = -RESISTI(3) / INDUCTI(3)
        ! --- LINHA 7 ---
        A(7,8)   =  1.0D0 / COMPLII(4)
        ! --- LINHA 8 ---
        A(8,3)   =  1.0D0 / INDUCTI(4)
        A(8,7)   = -1.0D0 / INDUCTI(4)
        A(8,8)   = -RESISTI(4) / INDUCTI(4)
        ! --- LINHA 9 ---
        A(9,10)  =  1.0D0 / COMPLII(5)
        A(9,12)  = -1.0D0 / COMPLII(5)
        A(9,14)  = -1.0D0 / COMPLII(5)
        ! --- LINHA 10 ---
        A(10,1)  = -1.0D0 / INDUCTI(5)
        A(10,9)  = -1.0D0 / INDUCTI(5)
        A(10,10) = -RESISTI(5) / INDUCTI(5)
        ! --- LINHA 11 ---
        A(11,12) =  1.0D0 / COMPLII(6)
        ! --- LINHA 12 ---
        A(12,9)  =  1.0D0 / INDUCTI(6)
        A(12,11) = -1.0D0 / INDUCTI(6)
        A(12,12) = -RESISTI(6) / INDUCTI(6)
        ! --- LINHA 13 ---
        A(13,14) =  1.0D0 / COMPLII(7)
        ! --- LINHA 14 ---
        A(14,9)  =  1.0D0 / INDUCTI(7)
        A(14,13) = -1.0D0 / INDUCTI(7)
        A(14,14) = -RESISTI(7) / INDUCTI(7)
        ! --- LINHA 15 ---
        A(15,16) =  1.0D0 / COMPLII(8)
        A(15,18) = -1.0D0 / COMPLII(8)
        A(15,20) = -1.0D0 / COMPLII(8)
        ! --- LINHA 16 ---
        A(16,1)  =  1.0D0 / INDUCTI(8)
        A(16,15) = -1.0D0 / INDUCTI(8)
        A(16,18) = -RESISTI(8) / INDUCTI(8)
        A(16,20) = -RESISTI(8) / INDUCTI(8)
        ! --- LINHA 17 ---
        A(17,18) =  1.0D0 / COMPLII(9)
        A(17,42) = -1.0D0 / COMPLII(9)
        A(17,44) = -1.0D0 / COMPLII(9)
        ! --- LINHA 18 ---
        A(18,15) =  1.0D0 / INDUCTI(9)
        A(18,17) = -1.0D0 / INDUCTI(9)
        A(18,42) = -RESISTI(9) / INDUCTI(9)
        A(18,44) = -RESISTI(9) / INDUCTI(9)
        ! --- LINHA 19 ---
        A(19,20) =  1.0D0 / COMPLII(10)
        A(19,22) = -1.0D0 / COMPLII(10)
        A(19,24) = -1.0D0 / COMPLII(10)
        ! --- LINHA 20 ---
        A(20,15) =  1.0D0 / INDUCTI(10)
        A(20,19) = -1.0D0 / INDUCTI(10)
        A(20,22) = -RESISTI(10) / INDUCTI(10)
        A(20,24) = -RESISTI(10) / INDUCTI(10)
        ! --- LINHA 21 ---
        A(21,22) =  1.0D0 / COMPLII(11)
        A(21,26) = -1.0D0 / COMPLII(11)
        A(21,28) = -1.0D0 / COMPLII(11)
        ! --- LINHA 22 ---
        A(22,19) =  1.0D0 / INDUCTI(11)
        A(22,21) = -1.0D0 / INDUCTI(11)
        A(22,26) = -RESISTI(11) / INDUCTI(11)
        A(22,28) = -RESISTI(11) / INDUCTI(11)
        ! --- LINHA 23 ---
        A(23,24) =  1.0D0 / COMPLII(12)
        A(23,38) = -1.0D0 / COMPLII(12)
        A(23,40) = -1.0D0 / COMPLII(12)
        ! --- LINHA 24 ---
        A(24,19) =  1.0D0 / INDUCTI(12)
        A(24,23) = -1.0D0 / INDUCTI(12)
        A(24,38) = -RESISTI(12) / INDUCTI(12)
        A(24,40) = -RESISTI(12) / INDUCTI(12)
        ! --- LINHA 25 ---
        A(25,26) =  1.0D0 / COMPLII(13)
        ! --- LINHA 26 ---
        A(26,21) =  1.0D0 / INDUCTI(13)
        A(26,25) = -1.0D0 / INDUCTI(13)
        A(26,26) = -RESISTI(13) / INDUCTI(13)
        ! --- LINHA 27 ---
        A(27,28) =  1.0D0 / COMPLII(14)
        A(27,30) = -1.0D0 / COMPLII(14)
        A(27,32) = -1.0D0 / COMPLII(14)
        ! --- LINHA 28 ---
        A(28,21) =  1.0D0 / INDUCTI(14)
        A(28,27) = -1.0D0 / INDUCTI(14)
        A(28,30) = -RESISTI(14) / INDUCTI(14)
        A(28,32) = -RESISTI(14) / INDUCTI(14)
        ! --- LINHA 29 ---
        A(29,30) =  1.0D0 / COMPLII(15)
        ! --- LINHA 30 ---
        A(30,27) =  1.0D0 / INDUCTI(15)
        A(30,29) = -1.0D0 / INDUCTI(15)
        A(30,30) = -RESISTI(15) / INDUCTI(15)
        ! --- LINHA 31 ---
        A(31,32) =  1.0D0 / COMPLII(16)
        A(31,34) = -1.0D0 / COMPLII(16)
        A(31,36) = -1.0D0 / COMPLII(16)
        ! --- LINHA 32 ---
        A(32,27) =  1.0D0 / INDUCTI(16)
        A(32,31) = -1.0D0 / INDUCTI(16)
        A(32,34) = -RESISTI(16) / INDUCTI(16)
        A(32,36) = -RESISTI(16) / INDUCTI(16)
        ! --- LINHA 33 ---
        A(33,34) =  1.0D0 / COMPLII(17)
        ! --- LINHA 34 ---
        A(34,31) =  1.0D0 / INDUCTI(17)
        A(34,33) = -1.0D0 / INDUCTI(17)
        A(34,34) = -RESISTI(17) / INDUCTI(17)
        ! --- LINHA 35 ---
        A(35,36) =  1.0D0 / COMPLII(18)
        ! --- LINHA 36 ---
        A(36,31) =  1.0D0 / INDUCTI(18)
        A(36,35) = -1.0D0 / INDUCTI(18)
        A(36,36) = -RESISTI(18) / INDUCTI(18)
        ! --- LINHA 37 ---
        A(37,38) =  1.0D0 / COMPLII(19)
        ! --- LINHA 38 ---
        A(38,23) =  1.0D0 / INDUCTI(19)
        A(38,37) = -1.0D0 / INDUCTI(19)
        A(38,38) = -RESISTI(19) / INDUCTI(19)
        ! --- LINHA 39 ---
        A(39,40) =  1.0D0 / COMPLII(20)
        ! --- LINHA 40 ---
        A(40,23) =  1.0D0 / INDUCTI(20)
        A(40,39) = -1.0D0 / INDUCTI(20)
        A(40,40) = -RESISTI(20) / INDUCTI(20)
        ! --- LINHA 41 ---
        A(41,42) =  1.0D0 / COMPLII(21)
        A(41,50) = -1.0D0 / COMPLII(21)
        A(41,52) = -1.0D0 / COMPLII(21)
        ! --- LINHA 42 ---
        A(42,17) =  1.0D0 / INDUCTI(21)
        A(42,41) = -1.0D0 / INDUCTI(21)
        A(42,50) = -RESISTI(21) / INDUCTI(21)
        A(42,52) = -RESISTI(21) / INDUCTI(21)
        ! --- LINHA 43 ---
        A(43,44) =  1.0D0 / COMPLII(22)
        A(43,46) = -1.0D0 / COMPLII(22)
        A(43,48) = -1.0D0 / COMPLII(22)
        ! --- LINHA 44 ---
        A(44,17) =  1.0D0 / INDUCTI(22)
        A(44,43) = -1.0D0 / INDUCTI(22)
        A(44,46) = -RESISTI(22) / INDUCTI(22)
        A(44,48) = -RESISTI(22) / INDUCTI(22)
        ! --- LINHA 45 ---
        A(45,46) =  1.0D0 / COMPLII(23)
        ! --- LINHA 46 ---
        A(46,43) =  1.0D0 / INDUCTI(23)
        A(46,45) = -1.0D0 / INDUCTI(23)
        A(46,46) = -RESISTI(23) / INDUCTI(23)
        ! --- LINHA 47 ---
        A(47,48) =  1.0D0 / COMPLII(24)
        ! --- LINHA 48 ---
        A(48,43) =  1.0D0 / INDUCTI(24)
        A(48,47) = -1.0D0 / INDUCTI(24)
        A(48,48) = -RESISTI(24) / INDUCTI(24)
        ! --- LINHA 49 ---
        A(49,50) =  1.0D0 / COMPLII(25)
        A(49,66) = -1.0D0 / COMPLII(25)
        A(49,68) = -1.0D0 / COMPLII(25)
        ! --- LINHA 50 ---
        A(50,41) =  1.0D0 / INDUCTI(25)
        A(50,49) = -1.0D0 / INDUCTI(25)
        A(50,66) = -RESISTI(25) / INDUCTI(25)
        A(50,68) = -RESISTI(25) / INDUCTI(25)
        ! --- LINHA 51 ---
        A(51,52) =  1.0D0 / COMPLII(26)
        A(51,54) = -1.0D0 / COMPLII(26)
        A(51,56) = -1.0D0 / COMPLII(26)
        ! --- LINHA 52 ---
        A(52,41) =  1.0D0 / INDUCTI(26)
        A(52,51) = -1.0D0 / INDUCTI(26)
        A(52,54) = -RESISTI(26) / INDUCTI(26)
        A(52,56) = -RESISTI(26) / INDUCTI(26)
        ! --- LINHA 53 ---
        A(53,54) =  1.0D0 / COMPLII(27)
        ! --- LINHA 54 ---
        A(54,51) =  1.0D0 / INDUCTI(27)
        A(54,53) = -1.0D0 / INDUCTI(27)
        A(54,54) = -RESISTI(27) / INDUCTI(27)
        ! --- LINHA 55 ---
        A(55,56) =  1.0D0 / COMPLII(28)
        A(55,58) = -1.0D0 / COMPLII(28)
        A(55,60) = -1.0D0 / COMPLII(28)
        ! --- LINHA 56 ---
        A(56,51) =  1.0D0 / INDUCTI(28)
        A(56,55) = -1.0D0 / INDUCTI(28)
        A(56,58) = -RESISTI(28) / INDUCTI(28)
        A(56,60) = -RESISTI(28) / INDUCTI(28)
        ! --- LINHA 57 ---
        A(57,58) =  1.0D0 / COMPLII(29)
        ! --- LINHA 58 ---
        A(58,55) =  1.0D0 / INDUCTI(29)
        A(58,57) = -1.0D0 / INDUCTI(29)
        A(58,58) = -RESISTI(29) / INDUCTI(29)
        ! --- LINHA 59 ---
        A(59,60) =  1.0D0 / COMPLII(30)
        A(59,62) = -1.0D0 / COMPLII(30)
        A(59,64) = -1.0D0 / COMPLII(30)
        ! --- LINHA 60 ---
        A(60,55) =  1.0D0 / INDUCTI(30)
        A(60,59) = -1.0D0 / INDUCTI(30)
        A(60,62) = -RESISTI(30) / INDUCTI(30)
        A(60,64) = -RESISTI(30) / INDUCTI(30)
        ! --- LINHA 61 ---
        A(61,62) =  1.0D0 / COMPLII(31)
        ! --- LINHA 62 ---
        A(62,59) =  1.0D0 / INDUCTI(31)
        A(62,61) = -1.0D0 / INDUCTI(31)
        A(62,62) = -RESISTI(31) / INDUCTI(31)
        ! --- LINHA 63 ---
        A(63,64) =  1.0D0 / COMPLII(32)
        ! --- LINHA 64 ---
        A(64,59) =  1.0D0 / INDUCTI(32)
        A(64,63) = -1.0D0 / INDUCTI(32)
        A(64,64) = -RESISTI(32) / INDUCTI(32)
        ! --- LINHA 65 ---
        A(65,66) =  1.0D0 / COMPLII(33)
        ! --- LINHA 66 ---
        A(66,49) =  1.0D0 / INDUCTI(33)
        A(66,65) = -1.0D0 / INDUCTI(33)
        A(66,66) = -RESISTI(33) / INDUCTI(33)
        ! --- LINHA 67 ---
        A(67,68) =  1.0D0 / COMPLII(34)
        A(67,70) = -1.0D0 / COMPLII(34)
        A(67,72) = -1.0D0 / COMPLII(34)
        ! --- LINHA 68 ---
        A(68,49) =  1.0D0 / INDUCTI(34)
        A(68,67) = -1.0D0 / INDUCTI(34)
        A(68,70) = -RESISTI(34) / INDUCTI(34)
        A(68,72) = -RESISTI(34) / INDUCTI(34)
        ! --- LINHA 69 ---
        A(69,70) =  1.0D0 / COMPLII(35)
        A(69,80) = -1.0D0 / COMPLII(35)
        A(69,82) = -1.0D0 / COMPLII(35)
        ! --- LINHA 70 ---
        A(70,67) =  1.0D0 / INDUCTI(35)
        A(70,69) = -1.0D0 / INDUCTI(35)
        A(70,80) = -RESISTI(35) / INDUCTI(35)
        A(70,82) = -RESISTI(35) / INDUCTI(35)
        ! --- LINHA 71 ---
        A(71,72) =  1.0D0 / COMPLII(36)
        A(71,74) = -1.0D0 / COMPLII(36)
        A(71,76) = -1.0D0 / COMPLII(36)
        A(71,78) = -1.0D0 / COMPLII(36)
        ! --- LINHA 72 ---
        A(72,67) =  1.0D0 / INDUCTI(36)
        A(72,71) = -1.0D0 / INDUCTI(36)
        A(72,74) = -RESISTI(36) / INDUCTI(36)
        A(72,76) = -RESISTI(36) / INDUCTI(36)
        A(72,78) = -RESISTI(36) / INDUCTI(36)
        ! --- LINHA 73 ---
        A(73,74) =  1.0D0 / COMPLII(37)
        ! --- LINHA 74 ---
        A(74,71) =  1.0D0 / INDUCTI(37)
        A(74,73) = -1.0D0 / INDUCTI(37)
        A(74,74) = -RESISTI(37) / INDUCTI(37)
        ! --- LINHA 75 ---
        A(75,76) =  1.0D0 / COMPLII(38)
        ! --- LINHA 76 ---
        A(76,71) =  1.0D0 / INDUCTI(38)
        A(76,75) = -1.0D0 / INDUCTI(38)
        A(76,76) = -RESISTI(38) / INDUCTI(38)
        ! --- LINHA 77 ---
        A(77,78) =  1.0D0 / COMPLII(39)
        ! --- LINHA 78 ---
        A(78,71) =  1.0D0 / INDUCTI(39)
        A(78,77) = -1.0D0 / INDUCTI(39)
        A(78,78) = -RESISTI(39) / INDUCTI(39)
        ! --- LINHA 79 ---
        A(79,80) =  1.0D0 / COMPLII(40)
        ! --- LINHA 80 ---
        A(80,69) =  1.0D0 / INDUCTI(40)
        A(80,79) = -1.0D0 / INDUCTI(40)
        A(80,80) = -RESISTI(40) / INDUCTI(40)
        ! --- LINHA 81 ---
        A(81,82) =  1.0D0 / COMPLII(41)
        A(81,84) = -1.0D0 / COMPLII(41)
        A(81,86) = -1.0D0 / COMPLII(41)
        ! --- LINHA 82 ---
        A(82,69) =  1.0D0 / INDUCTI(41)
        A(82,81) = -1.0D0 / INDUCTI(41)
        A(82,84) = -RESISTI(41) / INDUCTI(41)
        A(82,86) = -RESISTI(41) / INDUCTI(41)
        ! --- LINHA 83 ---
        A(83,84) =  1.0D0 / COMPLII(42)
        ! --- LINHA 84 ---
        A(84,81) =  1.0D0 / INDUCTI(42)
        A(84,83) = -1.0D0 / INDUCTI(42)
        A(84,84) = -RESISTI(42) / INDUCTI(42)
        ! --- LINHA 85 ---
        A(85,86) =  1.0D0 / COMPLII(43)
        A(85,88) = -1.0D0 / COMPLII(43)
        A(85,90) = -1.0D0 / COMPLII(43)
        ! --- LINHA 86 ---
        A(86,81) =  1.0D0 / INDUCTI(43)
        A(86,85) = -1.0D0 / INDUCTI(43)
        A(86,88) = -RESISTI(43) / INDUCTI(43)
        A(86,90) = -RESISTI(43) / INDUCTI(43)
        ! --- LINHA 87 ---
        A(87,88) =  1.0D0 / COMPLII(44)
        ! --- LINHA 88 ---
        A(88,85) =  1.0D0 / INDUCTI(44)
        A(88,87) = -1.0D0 / INDUCTI(44)
        A(88,88) = -RESISTI(44) / INDUCTI(44)
        ! --- LINHA 89 ---
        A(89,90) =  1.0D0 / COMPLII(45)
        A(89,92) = -1.0D0 / COMPLII(45)
        A(89,94) = -1.0D0 / COMPLII(45)
        ! --- LINHA 90 ---
        A(90,85) =  1.0D0 / INDUCTI(45)
        A(90,89) = -1.0D0 / INDUCTI(45)
        A(90,92) = -RESISTI(45) / INDUCTI(45)
        A(90,94) = -RESISTI(45) / INDUCTI(45)
        ! --- LINHA 91 ---
        A(91,92) =  1.0D0 / COMPLII(46)
        ! --- LINHA 92 ---
        A(92,89) =  1.0D0 / INDUCTI(46)
        A(92,91) = -1.0D0 / INDUCTI(46)
        A(92,92) = -RESISTI(46) / INDUCTI(46)
        ! --- LINHA 93 ---
        A(93,94) =  1.0D0 / COMPLII(47)
        A(93,96) = -1.0D0 / COMPLII(47)
        A(93,98) = -1.0D0 / COMPLII(47)
        ! --- LINHA 94 ---
        A(94,89) =  1.0D0 / INDUCTI(47)
        A(94,93) = -1.0D0 / INDUCTI(47)
        A(94,96) = -RESISTI(47) / INDUCTI(47)
        A(94,98) = -RESISTI(47) / INDUCTI(47)
        ! --- LINHA 95 ---
        A(95,96) =  1.0D0 / COMPLII(48)
        A(95,100) = -1.0D0 / COMPLII(48)
        A(95,102) = -1.0D0 / COMPLII(48)
        ! --- LINHA 96 ---
        A(96,93) =  1.0D0 / INDUCTI(48)
        A(96,95) = -1.0D0 / INDUCTI(48)
        A(96,100) = -RESISTI(48) / INDUCTI(48)
        A(96,102) = -RESISTI(48) / INDUCTI(48)
        ! --- LINHA 97 ---
        A(97,98) =  1.0D0 / COMPLII(49)
        A(97,112) = -1.0D0 / COMPLII(49)
        A(97,114) = -1.0D0 / COMPLII(49)
        ! --- LINHA 98 ---
        A(98,93) =  1.0D0 / INDUCTI(49)
        A(98,97) = -1.0D0 / INDUCTI(49)
        A(98,112) = -RESISTI(49) / INDUCTI(49)
        A(98,114) = -RESISTI(49) / INDUCTI(49)
        ! --- LINHA 99 ---
        A(99,100) =  1.0D0 / COMPLII(50)
        A(99,104) = -1.0D0 / COMPLII(50)
        A(99,106) = -1.0D0 / COMPLII(50)
        ! --- LINHA 100 ---
        A(100,95) =  1.0D0 / INDUCTI(50)
        A(100,99) = -1.0D0 / INDUCTI(50)
        A(100,104) = -RESISTI(50) / INDUCTI(50)
        A(100,106) = -RESISTI(50) / INDUCTI(50)
        ! --- LINHA 101 ---
        A(101,102) =  1.0D0 / COMPLII(51)
        ! --- LINHA 102 ---
        A(102,95) =  1.0D0 / INDUCTI(51)
        A(102,101) = -1.0D0 / INDUCTI(51)
        A(102,102) = -RESISTI(51) / INDUCTI(51)
        ! --- LINHA 103 ---
        A(103,104) =  1.0D0 / COMPLII(52)
        A(103,108) = -1.0D0 / COMPLII(52)
        A(103,110) = -1.0D0 / COMPLII(52)
        ! --- LINHA 104 ---
        A(104,100) =  1.0D0 / INDUCTI(52)
        A(104,104) = -1.0D0 / INDUCTI(52)
        A(104,108) = -RESISTI(52) / INDUCTI(52)
        A(104,110) = -RESISTI(52) / INDUCTI(52)
        ! --- LINHA 105 ---
        A(105,106) =  1.0D0 / COMPLII(53)
        ! --- LINHA 106 ---
        A(106,99) =  1.0D0 / INDUCTI(53)
        A(106,105) = -1.0D0 / INDUCTI(53)
        A(106,106) = -RESISTI(53) / INDUCTI(53)
        ! --- LINHA 107 ---
        A(107,108) =  1.0D0 / COMPLII(54)
        ! --- LINHA 108 ---
        A(108,103) =  1.0D0 / INDUCTI(54)
        A(108,107) = -1.0D0 / INDUCTI(54)
        A(108,108) = -RESISTI(54) / INDUCTI(54)
        ! --- LINHA 109 ---
        A(109,110) =  1.0D0 / COMPLII(55)
        ! --- LINHA 110 ---
        A(110,103) =  1.0D0 / INDUCTI(55)
        A(110,109) = -1.0D0 / INDUCTI(55)
        A(110,110) = -RESISTI(55) / INDUCTI(55)
        ! --- LINHA 111 ---
        A(111,112) =  1.0D0 / COMPLII(56)
        A(111,116) = -1.0D0 / COMPLII(56)
        A(111,118) = -1.0D0 / COMPLII(56)
        ! --- LINHA 112 ---
        A(112,97) =  1.0D0 / INDUCTI(56)
        A(112,111) = -1.0D0 / INDUCTI(56)
        A(112,116) = -RESISTI(56) / INDUCTI(56)
        A(112,118) = -RESISTI(56) / INDUCTI(56)
        ! --- LINHA 113 ---
        A(113,114) =  1.0D0 / COMPLII(57)
        ! --- LINHA 114 ---
        A(114,97) =  1.0D0 / INDUCTI(57)
        A(114,113) = -1.0D0 / INDUCTI(57)
        A(114,114) = -RESISTI(57) / INDUCTI(57)
        ! --- LINHA 115 ---
        A(115,116) =  1.0D0 / COMPLII(58)
        A(115,120) = -1.0D0 / COMPLII(58)
        A(115,122) = -1.0D0 / COMPLII(58)
        ! --- LINHA 116 ---
        A(116,111) =  1.0D0 / INDUCTI(58)
        A(116,115) = -1.0D0 / INDUCTI(58)
        A(116,120) = -RESISTI(58) / INDUCTI(58)
        A(116,122) = -RESISTI(58) / INDUCTI(58)
        ! --- LINHA 117 ---
        A(117,118) =  1.0D0 / COMPLII(59)
        ! --- LINHA 118 ---
        A(118,111) =  1.0D0 / INDUCTI(59)
        A(118,117) = -1.0D0 / INDUCTI(59)
        A(118,118) = -RESISTI(59) / INDUCTI(59)
        ! --- LINHA 119 ---
        A(119,120) =  1.0D0 / COMPLII(60)
        ! --- LINHA 120 ---
        A(120,115) =  1.0D0 / INDUCTI(60)
        A(120,119) = -1.0D0 / INDUCTI(60)
        A(120,120) = -RESISTI(60) / INDUCTI(60)
        ! --- LINHA 121 ---
        A(121,122) =  1.0D0 / COMPLII(61)
        ! --- LINHA 122 ---
        A(122,115) =  1.0D0 / INDUCTI(61)
        A(122,121) = -1.0D0 / INDUCTI(61)
        A(122,122) = -RESISTI(61) / INDUCTI(61)
        
        ! Vetor b: condições de contorno usando valores recebidos
        b(1) = PIN 
        b(2) = QIN 
    END SUBROUTINE FORMAB

    ! --------------------------------------------------------------------------
    ! MIVEULER: Integração pelo método de Euler explícito
    ! x_new = x_old + HT * (A @ x_old + b)
    ! --------------------------------------------------------------------------
    SUBROUTINE MIVEULER(dim, A, b, HT, IT, x_old, x_new)
        INTEGER, INTENT(IN) :: IT
        INTEGER, INTENT(IN)  :: dim
        REAL(8), INTENT(IN)  :: A(dim, dim), b(dim), HT, x_old(dim)
        REAL(8), INTENT(OUT) :: x_new(dim)

        REAL(8) :: dx(dim)
        
        ! dx = A * x_old + b  (multiplicação matriz-vetor)
        dx  = MATMUL(A, x_old) + b
        x_new = x_old + HT * dx
    END SUBROUTINE MIVEULER

END MODULE routines_mod


! ==============================================================================
! PROGRAMA PRINCIPAL
! ==============================================================================
PROGRAM ZERODMOD
    USE routines_mod
    IMPLICIT NONE

    ! Parâmetros da simulação
    INTEGER, PARAMETER :: N_FIXED = 61
    INTEGER, PARAMETER :: DIM     = 2 * N_FIXED
    INTEGER, PARAMETER :: N_MON   = 10

    REAL(8), PARAMETER :: DINA_TO_MMHG = 1D0/1333.22D0

    ! Variáveis
    INTEGER :: n, IT, NT, k, j, iunit, ios
    REAL(8) :: HT, TT, MI, RO, PI

    REAL(8), ALLOCATABLE :: COMPLI(:), RAIOI(:), ESPHI(:), MEY(:)
    REAL(8) :: RESISTI(N_FIXED), INDUCTI(N_FIXED), COMPLII(N_FIXED)

    REAL(8) :: A(DIM, DIM), b(DIM)
    REAL(8) :: x_old(DIM), x_new(DIM)

    REAL(8) :: PIN, QIN

    ! Segmentos monitorados (numeração 1-based, igual ao Python original)
    INTEGER :: indices_mon(N_MON)
    DATA indices_mon / 16, 19, 20, 23, 24, 30, 48, 49, 52, 58 /

    ! Array de resultados: (NT+1) linhas x (3 + N_MON*2) colunas
    INTEGER, PARAMETER :: NT_MAX = 100001  ! 3 segundos com HT=0.0001
    REAL(8) :: resultados(NT_MAX, 3 + N_MON*2)

    CHARACTER(LEN=20) :: hdr
    INTEGER :: col

    WRITE(*,*) 'Iniciando simulacao...'

    ! ------------------------------------------------------------------
    ! 1. Ler dados
    ! ------------------------------------------------------------------
    CALL DADOSMOD('DADOSMOD.txt', n, COMPLI, RAIOI, ESPHI, MEY, MI, RO, PI)
    WRITE(*,'(A,I4,A)') 'Trechos lidos: ', n

    IF (n /= N_FIXED) THEN
        WRITE(*,'(A,I4,A)') 'AVISO: codigo FORMAB tem indices fixos para n=61. Lido n=', n, '.'
    END IF

    ! ------------------------------------------------------------------
    ! 2. Calcular geometria / propriedades
    ! ------------------------------------------------------------------
    CALL GEOMAT(n, COMPLI, RAIOI, ESPHI, MEY, MI, RO, PI, RESISTI, INDUCTI, COMPLII)
    WRITE(*,*) 'Dados Geometricos Calculados.'

    ! ------------------------------------------------------------------
    ! 3. Configuração temporal
    ! ------------------------------------------------------------------
    HT = 0.0001
    TT = 10.0D0  ! Simular 3 segundos (3 ciclos cardíacos)
    NT = INT(TT / HT)   
    WRITE(*,'(A,I6,A,F5.1,A)') 'Integracao temporal (', NT, ' passos) para ', TT, ' segundos...'

    IF (NT + 1 > NT_MAX) THEN
        WRITE(*,*) 'ERRO: NT+1 excede NT_MAX. Ajuste o parametro NT_MAX no codigo.'
        STOP
    END IF

    ! ------------------------------------------------------------------
    ! 4. Inicialização
    ! ------------------------------------------------------------------
    A(:,:)    = 0.0D0
    b(:)      = 0.0D0
    x_old(:)  = 0.0D0
    x_new(:)  = 0.0D0
    resultados(:,:) = 0.0D0

    ! Preencher coluna de tempo (t = 0, HT, 2*HT, ..., TT)
    DO IT = 0, NT
        resultados(IT + 1, 1) = DBLE(IT) * HT
    END DO

    ! Estado inicial (t=0): calcular PIN/QIN
    CALL PULSOPQ(0, HT, PI, PIN, QIN)
    
    resultados(1, 2) = PIN
    resultados(1, 3) = QIN

    ! Preencher pontos monitorados para t=0
    DO k = 1, N_MON
        j = indices_mon(k)
        resultados(1, 3 + (k-1)*2 + 1) = x_old(2*j - 1)   ! Pressão
        resultados(1, 3 + (k-1)*2 + 2) = x_old(2*j)        ! Vazão
    END DO

    ! ------------------------------------------------------------------
    ! 5. Loop principal de simulação
    ! ------------------------------------------------------------------
    DO IT = 0, NT - 1
        ! Calcular PIN/QIN para tempo atual
        CALL PULSOPQ(IT, HT, PI, PIN, QIN)
        
        ! Armazenar PIN/QIN para tempo atual (exceto t=0 já feito)
        IF (IT > 0) THEN
            resultados(IT+1, 2) = PIN
            resultados(IT+1, 3) = QIN
        END IF

        ! Montar sistema com PIN/QIN calculados
        CALL FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PIN, QIN, A, b)
        CALL MIVEULER(DIM, A, b, HT, IT, x_old, x_new)

        ! Armazenar pontos monitorados para próximo tempo
        DO k = 1, N_MON
            j = indices_mon(k)
            resultados(IT + 2, 3 + (k-1)*2 + 1) = x_new(2*j - 1) 
            resultados(IT + 2, 3 + (k-1)*2 + 2) = x_new(2*j)
        END DO

        x_old(:) = x_new(:)
    END DO

    ! Armazenar último passo (t=TT)
    CALL PULSOPQ(NT, HT, PI, PIN, QIN)
    resultados(NT+1, 2) = PIN
    resultados(NT+1, 3) = QIN

    WRITE(*,*) 'Simulacao concluida.'

    ! ------------------------------------------------------------------
    ! 6. Salvar resultados em CSV
    ! ------------------------------------------------------------------
    iunit = 30
    OPEN(UNIT=iunit, FILE='resultados_simulacao.csv', STATUS='REPLACE', &
         ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
        WRITE(*,*) 'ERRO: Nao foi possivel criar resultados_simulacao.csv'
        STOP
    END IF

    ! Cabeçalho atualizado
    WRITE(iunit, '(A)', ADVANCE='NO') 'tempo,PIN,QIN'
    DO k = 1, N_MON
        j = indices_mon(k)
        WRITE(iunit, '(A,I0,A,I0)', ADVANCE='NO') ',P_', j, ',Q_', j
    END DO
    WRITE(iunit, *)   ! newline

    ! Dados
    DO IT = 1, NT + 1
        WRITE(iunit, '(ES15.6E2)', ADVANCE='NO') resultados(IT, 1)
        DO col = 2, 3 + N_MON*2
            WRITE(iunit, '(A,ES15.6E2)', ADVANCE='NO') ',', resultados(IT, col)
        END DO
        WRITE(iunit, *)
    END DO

    CLOSE(iunit)
    WRITE(*,*) 'Resultados salvos: resultados_simulacao.csv'

    ! ------------------------------------------------------------------
    ! 7. Salvar informações da simulação
    ! ------------------------------------------------------------------
    iunit = 31
    OPEN(UNIT=iunit, FILE='info_simulacao.csv', STATUS='REPLACE', ACTION='WRITE')
    WRITE(iunit,'(A)')       'Parametro,Valor'
    WRITE(iunit,'(A,ES15.6E2)') 'HT,', HT
    WRITE(iunit,'(A,ES15.6E2)') 'TT,', TT
    WRITE(iunit,'(A,I6)')    'NT,', NT
    WRITE(iunit,'(A,I6)')    'n,', n
    WRITE(iunit,'(A,ES15.6E2)') 'MI,', MI
    WRITE(iunit,'(A,ES15.6E2)') 'RO,', RO
    WRITE(iunit,'(A,ES15.6E2)') 'PI,', PI
    WRITE(iunit,'(A)', ADVANCE='NO') 'Indices_monitorados,'
    DO k = 1, N_MON
        IF (k < N_MON) THEN
            WRITE(iunit,'(I0,A)', ADVANCE='NO') indices_mon(k), ','
        ELSE
            WRITE(iunit,'(I0)') indices_mon(k)
        END IF
    END DO
    CLOSE(iunit)
    WRITE(*,*) 'Info salva: info_simulacao.csv'

    ! Liberar memória alocada
    IF (ALLOCATED(COMPLI)) DEALLOCATE(COMPLI)
    IF (ALLOCATED(RAIOI))  DEALLOCATE(RAIOI)
    IF (ALLOCATED(ESPHI))  DEALLOCATE(ESPHI)
    IF (ALLOCATED(MEY))    DEALLOCATE(MEY)

END PROGRAM ZERODMOD
