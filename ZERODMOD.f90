module routines_p
    implicit none
    public :: DADOSMOD, GEOMAT, PULSOPQ, FORMAB, MIVEULER
    private
    
contains

    subroutine DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO, PI
        real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
        integer :: i, ios
        real(kind=8) :: t_MI, t_RO, t_PI

        open(unit=10, file='DADOSMOD.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop "ERRO: DADOSMOD.txt nao encontrado."
        
        read(10,*) n
        
        allocate(COMPI(n), RAIOI(n), ESPHI(n), MEY(n))
        do i = 1, n
            read(10,*) COMPI(i), RAIOI(i), ESPHI(i), MEY(i), t_MI, t_RO, t_PI
            if (i == 1) then
                MI = t_MI
                RO = t_RO
                PI = t_PI
            endif
        end do
        close(10)
    end subroutine DADOSMOD

    subroutine GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI, RESISTI, INDUCTI, COMPLII)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO, PI
        real(kind=8) :: COMPI(n), RAIOI(n), ESPHI(n), MEY(n)
        real(kind=8) :: RESISTI(n), INDUCTI(n), COMPLII(n)
        real(kind=8) :: raio_quad, raio_cub, raio_quarta
        integer :: m, ios

        do m = 1, n
            raio_quad = RAIOI(m)**2
            raio_cub = RAIOI(m)**3
            raio_quarta = RAIOI(m)**4
            
            RESISTI(m) = (8.0d0 * MI * COMPI(m)) / (PI * raio_quarta)
            INDUCTI(m) = (9.0d0 * RO * COMPI(m)) / (4.0d0 * PI * raio_quad)
            COMPLII(m) = (3.0d0 * PI * COMPI(m) * raio_cub) / (2.0d0 * MEY(m) * ESPHI(m))
 
        end do

        open(unit=20, file='GEOMAT_R.txt', status='replace', action='write', iostat=ios)
        
        write(20,'(A,F12.6)') 'MI: ', MI
        write(20,'(A,F12.6)') 'RO: ', RO
        write(20,'(A,F12.8)') 'Pi: ', PI
        write(20,'(A,I4)')    'N : ', n
        write(20,'(A)') '================================================================'
        write(20,'(A)') ''
        write(20,'(A8,5A15)') 'Trecho', 'COMPI', 'RAIOI', 'ESPHI', 'MEY', 'RESISTI'
        write(20,'(A8,5A15)') '', '', '', '', '', 'INDUCTI'
        write(20,'(A8,5A15)') '', '', '', '', '', 'COMPLII'
        write(20,'(A)') '----------------------------------------------------------------'
        
        do m = 1, n
            write(20,'(I8,5ES15.6)') m, COMPI(m), RAIOI(m), ESPHI(m), MEY(m), RESISTI(m)
            write(20,'(8X,4(15X),ES15.6)') INDUCTI(m)
            write(20,'(8X,4(15X),ES15.6)') COMPLII(m)
            write(20,'(A)') ''
        end do
        
        write(20,'(A)') '================================================================'
        close(20)
    end subroutine GEOMAT

    subroutine PULSOPQ(IT, HT, PI, PIN, QIN)
        implicit none
        integer :: IT
        real(kind=8) :: HT, PI
        real(kind=8) :: PIN, QIN
        real(kind=8) :: TS, TCC, PS, PD, QMAX, TI, A, B

        TS = 0.3d0
        TCC = 1.0d0
        PS = 120.0d0
        PD = 80.0d0
        QMAX = 300.0d0
        
        TI = IT * HT
        
        if (TI <= TS) then
            PIN = (PD + ((PS - PD) / TS) * TI)* 1333.22
            QIN = QMAX * (SIN(PI * TI / TS))**2
        endif

        if (TI > TS .and. TI <= TCC) then
            A = (PS - PD) / (TS - TCC)
            B = (PD * TS - PS * TCC) / (TS - TCC)
            PIN = (A * TI + B) * 1333.22
            QIN = 0.0d0
        endif
    end subroutine PULSOPQ

    subroutine FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PI, IT, A, b)
        implicit none
        integer:: n, IT
        real(kind=8) :: HT, PI
        real(kind=8) :: RESISTI(n), INDUCTI(n), COMPLII(n)
        real(kind=8) :: A(122,122)
        real(kind=8) :: b(122)
        real(kind=8) :: PIN, QIN

        ! Matriz A
        A = 0.0d0
        ! --- LINHA 1 ---
        A(1,4)  =  1.0d0/COMPLII(1)
        A(1,10) =  1.0d0/COMPLII(1)
        A(1,16) =  1.0d0/COMPLII(1)
        
        ! --- LINHA 2 ---
        A(2,1)  =  1.0d0/INDUCTI(1)
        A(2,4)  = -RESISTI(1)/INDUCTI(1)
        A(2,10) = -RESISTI(1)/INDUCTI(1)
        A(2,16) = -RESISTI(1)/INDUCTI(1)
        
        ! --- LINHA 3 ---
        A(3,4) =  1.0d0/COMPLII(2)
        A(3,6) = -1.0d0/COMPLII(2)
        A(3,8) = -1.0d0/COMPLII(2)
        
        ! --- LINHA 4 ---
        A(4,1) =  1.0d0/INDUCTI(2)
        A(4,3) = -1.0d0/INDUCTI(2)
        A(4,4) = -RESISTI(2)/INDUCTI(2)
        
        ! --- LINHA 5 ---
        A(5,6) = 1.0d0/COMPLII(3)
        
        ! --- LINHA 6 ---
        A(6,5) =  1.0d0/INDUCTI(3)
        A(6,6) = -RESISTI(3)/INDUCTI(3)
        
        ! --- LINHA 7 ---
        A(7,8) = 1.0d0/COMPLII(4)
        
        ! --- LINHA 8 ---
        A(8,7) =  1.0d0/INDUCTI(4)
        A(8,8) = -RESISTI(4)/INDUCTI(4)
        
        ! --- LINHA 9 ---
        A(9,10) =  1.0d0/COMPLII(5)
        A(9,12) = -1.0d0/COMPLII(5)
        A(9,14) = -1.0d0/COMPLII(5)
        
        ! --- LINHA 10 ---
        A(10,1)  = -1.0d0/INDUCTI(5)
        A(10,9)  = -1.0d0/INDUCTI(5)
        A(10,10) = -RESISTI(5)/INDUCTI(5)
        
        ! --- LINHA 11 ---
        A(11,12) = 1.0d0/COMPLII(6)
        
        ! --- LINHA 12 ---
        A(12,11) =  1.0d0/INDUCTI(6)
        A(12,12) = -RESISTI(6)/INDUCTI(6)
        
        ! --- LINHA 13 ---
        A(13,14) = 1.0d0/COMPLII(7)
        
        ! --- LINHA 14 ---
        A(14,13) =  1.0d0/INDUCTI(7)
        A(14,14) = -RESISTI(7)/INDUCTI(7)
        
        ! --- LINHA 15 ---
        A(15,16) =  1.0d0/COMPLII(8)
        A(15,18) = -1.0d0/COMPLII(8)
        A(15,20) = -1.0d0/COMPLII(8)
        
        ! --- LINHA 16 ---
        A(16,1)  =  1.0d0/INDUCTI(8)
        A(16,15) = -1.0d0/INDUCTI(8)
        A(16,18) = -RESISTI(8)/INDUCTI(8)
        A(16,20) = -RESISTI(8)/INDUCTI(8)
        
        ! --- LINHA 17 ---
        A(17,18) =  1.0d0/COMPLII(9)
        A(17,42) = -1.0d0/COMPLII(9)
        A(17,44) = -1.0d0/COMPLII(9)
        
        ! --- LINHA 18 ---
        A(18,15) =  1.0d0/INDUCTI(9)
        A(18,17) = -1.0d0/INDUCTI(9)
        A(18,42) = -RESISTI(9)/INDUCTI(9)
        A(18,44) = -RESISTI(9)/INDUCTI(9)
        
        ! --- LINHA 19 ---
        A(19,20) =  1.0d0/COMPLII(10)
        A(19,22) = -1.0d0/COMPLII(10)
        A(19,24) = -1.0d0/COMPLII(10)
        
        ! --- LINHA 20 ---
        A(20,15) =  1.0d0/INDUCTI(10)
        A(20,19) = -1.0d0/INDUCTI(10)
        A(20,22) = -RESISTI(10)/INDUCTI(10)
        A(20,24) = -RESISTI(10)/INDUCTI(10)
        
        ! --- LINHA 21 ---
        A(21,22) =  1.0d0/COMPLII(11)
        A(21,26) = -1.0d0/COMPLII(11)
        A(21,28) = -1.0d0/COMPLII(11)
        
        ! --- LINHA 22 ---
        A(22,19) =  1.0d0/INDUCTI(11)
        A(22,21) = -1.0d0/INDUCTI(11)
        A(22,26) = -RESISTI(11)/INDUCTI(11)
        A(22,28) = -RESISTI(11)/INDUCTI(11)
        
        ! --- LINHA 23 ---
        A(23,24) =  1.0d0/COMPLII(12)
        A(23,38) = -1.0d0/COMPLII(12)
        A(23,40) = -1.0d0/COMPLII(12)
        
        ! --- LINHA 24 ---
        A(24,19) =  1.0d0/INDUCTI(12)
        A(24,23) = -1.0d0/INDUCTI(12)
        A(24,38) = -RESISTI(12)/INDUCTI(12)
        A(24,40) = -RESISTI(12)/INDUCTI(12)
        
        ! --- LINHA 25 ---
        A(25,26) = 1.0d0/COMPLII(13)
        
        ! --- LINHA 26 ---
        A(26,25) =  1.0d0/INDUCTI(13)
        A(26,26) = -RESISTI(13)/INDUCTI(13)
        
        ! --- LINHA 27 ---
        A(27,28) =  1.0d0/COMPLII(14)
        A(27,30) = -1.0d0/COMPLII(14)
        A(27,32) = -1.0d0/COMPLII(14)
        
        ! --- LINHA 28 ---
        A(28,21) =  1.0d0/INDUCTI(14)
        A(28,27) = -1.0d0/INDUCTI(14)
        A(28,30) = -RESISTI(14)/INDUCTI(14)
        A(28,32) = -RESISTI(14)/INDUCTI(14)
        
        ! --- LINHA 29 ---
        A(29,30) = 1.0d0/COMPLII(15)
        
        ! --- LINHA 30 ---
        A(30,29) =  1.0d0/INDUCTI(15)
        A(30,30) = -RESISTI(15)/INDUCTI(15)
        
        ! --- LINHA 31 ---
        A(31,32) =  1.0d0/COMPLII(16)
        A(31,34) = -1.0d0/COMPLII(16)
        A(31,36) = -1.0d0/COMPLII(16)
        
        ! --- LINHA 32 ---
        A(32,27) =  1.0d0/INDUCTI(16)
        A(32,31) = -1.0d0/INDUCTI(16)
        A(32,34) = -RESISTI(16)/INDUCTI(16)
        A(32,36) = -RESISTI(16)/INDUCTI(16)
        
        ! --- LINHA 33 ---
        A(33,34) = 1.0d0/COMPLII(17)
        
        ! --- LINHA 34 ---
        A(34,33) =  1.0d0/INDUCTI(17)
        A(34,34) = -RESISTI(17)/INDUCTI(17)
        
        ! --- LINHA 35 ---
        A(35,36) = 1.0d0/COMPLII(18)
        
        ! --- LINHA 36 ---
        A(36,35) =  1.0d0/INDUCTI(18)
        A(36,36) = -RESISTI(18)/INDUCTI(18)
        
        ! --- LINHA 37 ---
        A(37,38) = 1.0d0/COMPLII(19)
        
        ! --- LINHA 38 ---
        A(38,37) =  1.0d0/INDUCTI(19)
        A(38,38) = -RESISTI(19)/INDUCTI(19)
        
        ! --- LINHA 39 ---
        A(39,40) = 1.0d0/COMPLII(20)
        
        ! --- LINHA 40 ---
        A(40,39) =  1.0d0/INDUCTI(20)
        A(40,40) = -RESISTI(20)/INDUCTI(20)
        
        ! --- LINHA 41 ---
        A(41,42) =  1.0d0/COMPLII(21)
        A(41,50) = -1.0d0/COMPLII(21)
        A(41,52) = -1.0d0/COMPLII(21)
        
        ! --- LINHA 42 ---
        A(42,17) =  1.0d0/INDUCTI(21)
        A(42,41) = -1.0d0/INDUCTI(21)
        A(42,50) = -RESISTI(21)/INDUCTI(21)
        A(42,52) = -RESISTI(21)/INDUCTI(21)
        
        ! --- LINHA 43 ---
        A(43,44) =  1.0d0/COMPLII(22)
        A(43,46) = -1.0d0/COMPLII(22)
        A(43,48) = -1.0d0/COMPLII(22)
        
        ! --- LINHA 44 ---
        A(44,17) =  1.0d0/INDUCTI(22)
        A(44,43) = -1.0d0/INDUCTI(22)
        A(44,46) = -RESISTI(22)/INDUCTI(22)
        A(44,48) = -RESISTI(22)/INDUCTI(22)
        
        ! --- LINHA 45 ---
        A(45,46) = 1.0d0/COMPLII(23)
        
        ! --- LINHA 46 ---
        A(46,45) =  1.0d0/INDUCTI(23)
        A(46,46) = -RESISTI(23)/INDUCTI(23)
        
        ! --- LINHA 47 ---
        A(47,48) = 1.0d0/COMPLII(24)
        
        ! --- LINHA 48 ---
        A(48,47) =  1.0d0/INDUCTI(24)
        A(48,48) = -RESISTI(24)/INDUCTI(24)
        
        ! --- LINHA 49 ---
        A(49,50) =  1.0d0/COMPLII(25)
        A(49,66) = -1.0d0/COMPLII(25)
        A(49,68) = -1.0d0/COMPLII(25)
        
        ! --- LINHA 50 ---
        A(50,41) =  1.0d0/INDUCTI(25)
        A(50,49) = -1.0d0/INDUCTI(25)
        A(50,66) = -RESISTI(25)/INDUCTI(25)
        A(50,68) = -RESISTI(25)/INDUCTI(25)
        
        ! --- LINHA 51 ---
        A(51,52) =  1.0d0/COMPLII(26)
        A(51,54) = -1.0d0/COMPLII(26)
        A(51,56) = -1.0d0/COMPLII(26)
        
        ! --- LINHA 52 ---
        A(52,41) =  1.0d0/INDUCTI(26)
        A(52,51) = -1.0d0/INDUCTI(26)
        A(52,54) = -RESISTI(26)/INDUCTI(26)
        A(52,56) = -RESISTI(26)/INDUCTI(26)
        
        ! --- LINHA 53 ---
        A(53,54) = 1.0d0/COMPLII(27)
        
        ! --- LINHA 54 ---
        A(54,53) =  1.0d0/INDUCTI(27)
        A(54,54) = -RESISTI(27)/INDUCTI(27)
        
        ! --- LINHA 55 ---
        A(55,56) =  1.0d0/COMPLII(28)
        A(55,58) = -1.0d0/COMPLII(28)
        A(55,60) = -1.0d0/COMPLII(28)
        
        ! --- LINHA 56 ---
        A(56,51) =  1.0d0/INDUCTI(28)
        A(56,55) = -1.0d0/INDUCTI(28)
        A(56,58) = -RESISTI(28)/INDUCTI(28)
        A(56,60) = -RESISTI(28)/INDUCTI(28)
        
        ! --- LINHA 57 ---
        A(57,58) = 1.0d0/COMPLII(29)
        
        ! --- LINHA 58 ---
        A(58,57) =  1.0d0/INDUCTI(29)
        A(58,58) = -RESISTI(29)/INDUCTI(29)
        
        ! --- LINHA 59 ---
        A(59,60) =  1.0d0/COMPLII(30)
        A(59,62) = -1.0d0/COMPLII(30)
        A(59,64) = -1.0d0/COMPLII(30)
        
        ! --- LINHA 60 ---
        A(60,55) =  1.0d0/INDUCTI(30)
        A(60,59) = -1.0d0/INDUCTI(30)
        A(60,62) = -RESISTI(30)/INDUCTI(30)
        A(60,64) = -RESISTI(30)/INDUCTI(30)
        
        ! --- LINHA 61 ---
        A(61,62) = 1.0d0/COMPLII(31)
        
        ! --- LINHA 62 ---
        A(62,61) =  1.0d0/INDUCTI(31)
        A(62,62) = -RESISTI(31)/INDUCTI(31)
        
        ! --- LINHA 63 ---
        A(63,64) = 1.0d0/COMPLII(32)
        
        ! --- LINHA 64 ---
        A(64,63) =  1.0d0/INDUCTI(32)
        A(64,64) = -RESISTI(32)/INDUCTI(32)
        
        ! --- LINHA 65 ---
        A(65,66) = 1.0d0/COMPLII(33)
        
        ! --- LINHA 66 ---
        A(66,65) =  1.0d0/INDUCTI(33)
        A(66,66) = -RESISTI(33)/INDUCTI(33)
        
        ! --- LINHA 67 ---
        A(67,68) =  1.0d0/COMPLII(34)
        A(67,70) = -1.0d0/COMPLII(34)
        A(67,72) = -1.0d0/COMPLII(34)
        
        ! --- LINHA 68 ---
        A(68,49) =  1.0d0/INDUCTI(34)
        A(68,67) = -1.0d0/INDUCTI(34)
        A(68,70) = -RESISTI(34)/INDUCTI(34)
        A(68,72) = -RESISTI(34)/INDUCTI(34)
        
        ! --- LINHA 69 ---
        A(69,70) =  1.0d0/COMPLII(35)
        A(69,80) = -1.0d0/COMPLII(35)
        A(69,82) = -1.0d0/COMPLII(35)
        
        ! --- LINHA 70 ---
        A(70,67) =  1.0d0/INDUCTI(35)
        A(70,69) = -1.0d0/INDUCTI(35)
        A(70,80) = -RESISTI(35)/INDUCTI(35)
        A(70,82) = -RESISTI(35)/INDUCTI(35)
        
        ! --- LINHA 71 ---
        A(71,72) =  1.0d0/COMPLII(36)
        A(71,74) = -1.0d0/COMPLII(36)
        A(71,76) = -1.0d0/COMPLII(36)
        A(71,78) = -1.0d0/COMPLII(36)
        
        ! --- LINHA 72 ---
        A(72,67) =  1.0d0/INDUCTI(36)
        A(72,71) = -1.0d0/INDUCTI(36)
        A(72,74) = -RESISTI(36)/INDUCTI(36)
        A(72,76) = -RESISTI(36)/INDUCTI(36)
        A(72,78) = -RESISTI(36)/INDUCTI(36)
        
        ! --- LINHA 73 ---
        A(73,74) = 1.0d0/COMPLII(37)
        
        ! --- LINHA 74 ---
        A(74,73) =  1.0d0/INDUCTI(37)
        A(74,74) = -RESISTI(37)/INDUCTI(37)
        
        ! --- LINHA 75 ---
        A(75,76) = 1.0d0/COMPLII(38)
        
        ! --- LINHA 76 ---
        A(76,75) =  1.0d0/INDUCTI(38)
        A(76,76) = -RESISTI(38)/INDUCTI(38)
        
        ! --- LINHA 77 ---
        A(77,78) = 1.0d0/COMPLII(39)
        
        ! --- LINHA 78 ---
        A(78,77) =  1.0d0/INDUCTI(39)
        A(78,78) = -RESISTI(39)/INDUCTI(39)
        
        ! --- LINHA 79 ---
        A(79,80) = 1.0d0/COMPLII(40)
        
        ! --- LINHA 80 ---
        A(80,79) =  1.0d0/INDUCTI(40)
        A(80,80) = -RESISTI(40)/INDUCTI(40)
        
        ! --- LINHA 81 ---
        A(81,82) =  1.0d0/COMPLII(41)
        A(81,84) = -1.0d0/COMPLII(41)
        A(81,86) = -1.0d0/COMPLII(41)
        
        ! --- LINHA 82 ---
        A(82,69) =  1.0d0/INDUCTI(41)
        A(82,81) = -1.0d0/INDUCTI(41)
        A(82,84) = -RESISTI(41)/INDUCTI(41)
        A(82,86) = -RESISTI(41)/INDUCTI(41)
        
        ! --- LINHA 83 ---
        A(83,84) = 1.0d0/COMPLII(42)
        
        ! --- LINHA 84 ---
        A(84,83) =  1.0d0/INDUCTI(42)
        A(84,84) = -RESISTI(42)/INDUCTI(42)
        
        ! --- LINHA 85 ---
        A(85,86) =  1.0d0/COMPLII(43)
        A(85,88) = -1.0d0/COMPLII(43)
        A(85,90) = -1.0d0/COMPLII(43)
        
        ! --- LINHA 86 ---
        A(86,81) =  1.0d0/INDUCTI(43)
        A(86,85) = -1.0d0/INDUCTI(43)
        A(86,88) = -RESISTI(43)/INDUCTI(43)
        A(86,90) = -RESISTI(43)/INDUCTI(43)
        
        ! --- LINHA 87 ---
        A(87,88) = 1.0d0/COMPLII(44)
        
        ! --- LINHA 88 ---
        A(88,87) =  1.0d0/INDUCTI(44)
        A(88,88) = -RESISTI(44)/INDUCTI(44)
        
        ! --- LINHA 89 ---
        A(89,90) =  1.0d0/COMPLII(45)
        A(89,92) = -1.0d0/COMPLII(45)
        A(89,94) = -1.0d0/COMPLII(45)
        
        ! --- LINHA 90 ---
        A(90,85) =  1.0d0/INDUCTI(45)
        A(90,89) = -1.0d0/INDUCTI(45)
        A(90,92) = -RESISTI(45)/INDUCTI(45)
        A(90,94) = -RESISTI(45)/INDUCTI(45)
        
        ! --- LINHA 91 ---
        A(91,92) = 1.0d0/COMPLII(46) 
        
        ! --- LINHA 92 ---
        A(92,91) = 1.0d0/INDUCTI(46) 
        A(92,92) = -RESISTI(46)/INDUCTI(46)

        ! --- LINHA 93 ---
        A(93,94) = 1.0d0/COMPLII(47) 
        A(93,96) = -1.0d0/COMPLII(47) 
        A(93,98) = -1.0d0/COMPLII(47)

        ! --- LINHA 94 ---
        A(94,89) = 1.0d0/INDUCTI(47)
        A(94,93) = -1.0d0/INDUCTI(47)
        A(94,96) = -RESISTI(47)/INDUCTI(47)
        A(94,98) = -RESISTI(47)/INDUCTI(47)

        ! --- LINHA 95 ---
        A(95,96) = 1.0d0/COMPLII(48)
        A(95,100) = -1.0d0/COMPLII(48)
        A(95,102) = -1.0d0/COMPLII(48)

        ! --- LINHA 96 ---
        A(96,93) = 1.0d0/INDUCTI(48)
        A(96,95) = -1.0d0/INDUCTI(48)
        A(96,100) = -RESISTI(48)/INDUCTI(48)
        A(96,102) = -RESISTI(48)/INDUCTI(48)

        ! --- LINHA 97 ---
        A(97,98) = 1.0d0/COMPLII(49)
        A(97,112) = -1.0d0/COMPLII(49)
        A(97,114) = -1.0d0/COMPLII(49)

        ! --- LINHA 98 ---
        A(98,93) = 1.0d0/INDUCTI(49)
        A(98,97) = -1.0d0/INDUCTI(49)
        A(98,112) = -RESISTI(49)/INDUCTI(49)
        A(98,114) = -RESISTI(49)/INDUCTI(49)

        ! --- LINHA 99 ---
        A(99,100) = 1.0d0/COMPLII(50)
        A(99,104) = -1.0d0/COMPLII(50)
        A(99,106) = -1.0d0/COMPLII(50)

        ! --- LINHA 100 ---
        A(100,95) = 1.0d0/INDUCTI(50)
        A(100,99) = -1.0d0/INDUCTI(50)
        A(100,104) = -RESISTI(50)/INDUCTI(50)
        A(100,106) = -RESISTI(50)/INDUCTI(50)

        ! --- LINHA 101 ---
        A(101,102) = 1.0d0/COMPLII(51)
        
        ! --- LINHA 102 ---
        A(102,101) = 1.0d0/INDUCTI(51)
        A(102,102) = -RESISTI(51)/INDUCTI(51)

        ! --- LINHA 103 ---
        A(103,104) = 1.0d0/COMPLII(52)
        A(103,108) = -1.0d0/COMPLII(52)
        A(103,110) = -1.0d0/COMPLII(52)

        ! --- LINHA 104 ---
        A(104,99) = 1.0d0/INDUCTI(52)
        A(104,103) = -1.0d0/INDUCTI(52)
        A(104,108) = -RESISTI(52)/INDUCTI(52)
        A(104,110) = -RESISTI(52)/INDUCTI(52)

        ! --- LINHA 105 ---
        A(105,106) = 1.0d0/COMPLII(53)
        
        ! --- LINHA 106 ---
        A(106,105) = 1.0d0/INDUCTI(53)
        A(106,106) = -RESISTI(53)/INDUCTI(53)

        ! --- LINHA 107 ---
        A(107,108) = 1.0d0/COMPLII(54)
        
        ! --- LINHA 108 ---
        A(108,107) = 1.0d0/INDUCTI(54)
        A(108,108) = -RESISTI(54)/INDUCTI(54)

        ! --- LINHA 109 ---
        A(109,110) = 1.0d0/COMPLII(55)
        
        ! --- LINHA 110 ---
        A(110,109) = 1.0d0/INDUCTI(55)
        A(110,110) = -RESISTI(55)/INDUCTI(55)

        ! --- LINHA 111 ---
        A(111,112) = 1.0d0/COMPLII(56)
        A(111,116) = -1.0d0/COMPLII(56)
        A(111,118) = -1.0d0/COMPLII(56)

        ! --- LINHA 112 ---
        A(112,97) = 1.0d0/INDUCTI(56)
        A(112,111) = -1.0d0/INDUCTI(56)
        A(112,116) = -RESISTI(56)/INDUCTI(56)
        A(112,118) = -RESISTI(56)/INDUCTI(56)

        ! --- LINHA 113 ---
        A(113,114) = 1.0d0/COMPLII(57)
        
        ! --- LINHA 114 ---
        A(114,113) = 1.0d0/INDUCTI(57)
        A(114,114) = -RESISTI(57)/INDUCTI(57)

        ! --- LINHA 115 ---
        A(115,116) = 1.0d0/COMPLII(58)
        A(115,120) = -1.0d0/COMPLII(58)
        A(115,122) = -1.0d0/COMPLII(58)

        ! --- LINHA 116 ---
        A(116,111) = 1.0d0/INDUCTI(58)
        A(116,115) = -1.0d0/INDUCTI(58)
        A(116,120) = -RESISTI(58)/INDUCTI(58)
        A(116,122) = -RESISTI(58)/INDUCTI(58)

        ! --- LINHA 117 ---
        A(117,118) = 1.0d0/COMPLII(59)
        
        ! --- LINHA 118 ---
        A(118,117) = 1.0d0/INDUCTI(59)
        A(118,118) = -RESISTI(59)/INDUCTI(59)

        ! --- LINHA 119 ---
        A(119,120) = 1.0d0/COMPLII(60)
        
        ! --- LINHA 120 ---
        A(120,119) = 1.0d0/INDUCTI(60)
        A(120,120) = -RESISTI(60)/INDUCTI(60)

        ! --- LINHA 121 ---
        A(121,122) = 1.0d0/COMPLII(61)
        
        ! --- LINHA 122 ---
        A(122,121) = 1.0d0/INDUCTI(61)
        A(122,122) = -RESISTI(61)/INDUCTI(61)

        ! Vetor b
        b = 0.0d0
        call PULSOPQ(IT, HT, PI, PIN, QIN)
        b(1) = PIN
        b(2) = QIN
        print *, b

    end subroutine FORMAB

    ! Metodo de euler vetorial Ax + b
    subroutine MIVEULER(n, A, b, HT, x_old, x_new)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: HT
        real(kind=8), intent(in) :: A(122,122)
        real(kind=8), intent(in) :: b(122), x_old(122)
        real(kind=8), intent(out) :: x_new(122)
        real(kind=8) :: dx(122)
        integer :: i, j
        integer :: dim
        
        ! Dimensão efetiva do sistema
        dim = 2*n
        
        ! Verificação de segurança
        if (dim > 122) then
            print *, "ERRO: dimensão excede tamanho do array"
            stop
        endif
        
        ! Inicialização
        dx = 0.0d0
        
        ! Calcula dx = A*x_old + b
        do i = 1, dim
            do j = 1, dim
                dx(i) = dx(i) + A(i,j) * x_old(j)
            end do
            dx(i) = dx(i) + b(i)
        end do
        
        ! Método de Euler: x_new = x_old + HT * dx
        do i = 1, dim
            x_new(i) = x_old(i) + HT * dx(i)
        end do

    end subroutine MIVEULER

end module routines_p
    
program ZERODMOD
    use routines_p
    implicit none
    
    integer :: n, IT, NT, i, j, k, ios
    integer, dimension(4) :: indices
    real(kind=8) :: MI, RO, PI, HT, TT, tempo
    real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
    real(kind=8), dimension(:), allocatable :: RESISTI, INDUCTI, COMPLII
    real(kind=8), dimension(122,122) :: A
    real(kind=8), dimension(122) :: b, x_old, x_new
    real(kind=8), dimension(:,:), allocatable :: resultados
    
    ! Índices dos trechos que queremos monitorar
    indices = (/16, 30, 48, 49/)
    
    write(*,'(A)') '========================================================'
    write(*,'(A)') '           SIMULACAO ZERODMOD - INICIANDO'
    write(*,'(A)') '========================================================'
    write(*,*)
    
    call DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI)
    write(*,'(A,I4)') 'Numero de trechos lidos: ', n
    write(*,*)
    
    allocate(RESISTI(n), INDUCTI(n), COMPLII(n))
    call GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI, RESISTI, INDUCTI, COMPLII)
    write(*,'(A)') 'Propriedades geometricas calculadas e salvas.'
    write(*,*)
    
    HT = 0.01d0
    TT = 1.0d0
    NT = int(TT / HT)
    
    allocate(resultados(NT+1, 1+4*2))  ! tempo + 4 trechos × 2 variáveis (P e Q)
    
    x_old = 0.0d0
    x_new = 0.0d0
    resultados = 0.0d0
    
    write(*,'(A)') 'Iniciando simulacao temporal...'
    write(*,'(A,I8)') 'Numero de passos de tempo: ', NT
    write(*,*)
    
    do IT = 0, NT
        tempo = IT * HT
        
        ! Forma a matriz A e vetor b para este instante de tempo
        call FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PI, IT, A, b)
        
        ! Aplica método de Euler
        call MIVEULER(n, A, b, HT, x_old, x_new)
        
        ! Armazena resultados para os trechos de interesse
        resultados(IT+1, 1) = tempo
        do k = 1, 4
            j = indices(k)
            ! Pressão no trecho j está na posição 2*j-1
            ! Vazão no trecho j está na posição 2*j
            resultados(IT+1, 2 + (k-1)*2) = x_new(2*j - 1)  ! Pressão
            resultados(IT+1, 3 + (k-1)*2) = x_new(2*j)      ! Vazão
        end do
        
        ! Atualiza para próximo passo
        x_old = x_new
        
        if (mod(IT, max(1, NT/10)) == 0) then
            write(*,'(A,F8.2,A)') 'Progresso: ', (100.0d0*IT)/NT, '%'
        endif
    end do
    
    write(*,*)
    write(*,'(A)') 'Simulacao concluida!'
    write(*,*)
    
    ! Salva resultados
    open(unit=30, file='resultados_simulacao.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) stop "ERRO: Nao foi possivel criar arquivo de resultados."
    
    write(30,'(A)') 'tempo,P_16,Q_16,P_30,Q_30,P_48,Q_48,P_49,Q_49'
    
    do IT = 1, NT+1
        write(30,'(ES15.8,8(A,ES15.8))') resultados(IT,1), &
                                        (',', resultados(IT,i), i=2,9)
    end do
    
    close(30)
    
    write(*,'(A)') 'Resultados salvos em: resultados_simulacao.csv'
    write(*,*)
    
    ! Arquivo de informações
    open(unit=31, file='info_simulacao.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) stop "ERRO: Nao foi possivel criar arquivo de informacoes."
    
    write(31,'(A)') '================================================================'
    write(31,'(A)') '         INFORMACOES DA SIMULACAO ZERODMOD'
    write(31,'(A)') '================================================================'
    write(31,'(A,F12.8)') 'Passo de tempo (HT): ', HT
    write(31,'(A,F12.6)') 'Tempo total (TT)   : ', TT
    write(31,'(A,I10)')   'Num. de passos (NT): ', NT
    write(31,'(A,I6)')    'Num. de trechos (n): ', n
    write(31,'(A,F12.8)') 'Viscosidade (MI)   : ', MI
    write(31,'(A,F12.6)') 'Densidade (RO)     : ', RO
    write(31,'(A,F12.10)')'Pi                 : ', PI
    write(31,'(A)') 'Indices monitorados    : 16, 30, 48, 49'
    write(31,'(A)') '================================================================'
    write(31,'(A)') ''
    write(31,'(A)') 'Estrutura do vetor de estado:'
    write(31,'(A)') '  Posicao 2*j-1: Pressao no trecho j'
    write(31,'(A)') '  Posicao 2*j  : Vazao no trecho j'
    write(31,'(A)') ''
    write(31,'(A)') 'Colunas no arquivo CSV:'
    write(31,'(A)') '  1: tempo'
    write(31,'(A)') '  2-3:  P_16, Q_16'
    write(31,'(A)') '  4-5:  P_30, Q_30'
    write(31,'(A)') '  6-7:  P_48, Q_48'
    write(31,'(A)') '  8-9:  P_49, Q_49'
    write(31,'(A)') '================================================================'
    
    close(31)
    
    write(*,'(A)') 'Informacoes salvas em: info_simulacao.txt'
    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A)') '           SIMULACAO ZERODMOD - FINALIZADA'
    write(*,'(A)') '========================================================'
    
    deallocate(COMPI, RAIOI, ESPHI, MEY)
    deallocate(RESISTI, INDUCTI, COMPLII)
    deallocate(resultados)
    
end program ZERODMOD

