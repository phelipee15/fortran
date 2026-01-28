module routines_p
    implicit none
    public :: DADOSMOD, GEOMAT, FORMAB, MIVEULER, SALVAR_PROPRIEDADES, VISUALIZAR_MATRIZ
    private
    real(kind=8), parameter :: PI = 3.141592653589793d0
    
contains

    subroutine DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO
        real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
        integer :: i, ios
        real(kind=8) :: temp_MI, temp_RO  ! Variáveis temporárias para leitura

        open(unit=10, file='DADOSMOD.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel abrir o arquivo DADOSMOD.txt"
            stop
        endif

        ! Lê número de elementos
        read(10,*,iostat=ios) n
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel ler o numero de elementos"
            stop
        endif
        
        print *, "========================================="
        print *, "Numero de elementos: ", n
        print *, "========================================="

        allocate(COMPI(n), RAIOI(n), ESPHI(n), MEY(n))

        ! Lê todos os dados em um único loop (6 valores por linha)
        print *, "Lendo dados completos (6 valores por linha)..."
        do i = 1, n
            read(10,*,iostat=ios) COMPI(i), RAIOI(i), ESPHI(i), MEY(i), temp_MI, temp_RO
            if (ios /= 0) then
                print *, "ERRO: Falha ao ler linha ", i, " do arquivo"
                print *, "Verifique se o arquivo tem exatamente 6 valores por linha"
                stop
            endif
            
            ! Armazena MI e RO da PRIMEIRA linha (são constantes para todos os elementos)
            if (i == 1) then
                MI = temp_MI
                RO = temp_RO
            else
                ! Verificação de consistência (opcional mas recomendada)
                if (abs(temp_MI - MI) > 1.0d-12 .or. abs(temp_RO - RO) > 1.0d-6) then
                    print *, "AVISO: Valores de MI/RO inconsistentes na linha", i
                    print *, "  Esperado: MI =", MI, " RO =", RO
                    print *, "  Lido:     MI =", temp_MI, " RO =", temp_RO
                endif
            endif
            
            ! Impressão formatada dos 6 valores
            print '(A,I3,A,6F10.4)', "  Elemento ", i, ": ", &
                  COMPI(i), RAIOI(i), ESPHI(i), MEY(i), temp_MI, temp_RO
        end do

        close(10)

        print *, "========================================="
        print *, "Propriedades do Material (constantes):"
        print *, "  MI (permeabilidade): ", MI
        print *, "  RO (densidade):      ", RO
        print *, "========================================="
    end subroutine DADOSMOD

    subroutine GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, RESISTI, INDUCTI, COMPLII)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO
        real(kind=8), dimension(n) :: COMPI, RAIOI, ESPHI, MEY
        real(kind=8), dimension(n) :: RESISTI, INDUCTI, COMPLII
        integer :: m
        real(kind=8) :: raio_quad, raio_cub, raio_quarta
        
        print *, "Calculando propriedades dos elementos..."
        print *, ""
        
        do m = 1, n
            raio_quad = RAIOI(m)**2
            raio_cub = RAIOI(m)**3
            raio_quarta = RAIOI(m)**4
            
            ! Resistência (considerando efeito do número de espiras implícito em MEY)
            RESISTI(m) = (8.0d0 * MI * COMPI(m)) / (PI * raio_quarta)
            
            ! Indutância
            INDUCTI(m) = (9.0d0 * RO * COMPI(m)) / (4.0d0 * PI * raio_quad)
            
            ! Compliância (flexibilidade) - MEY agora vem diretamente do arquivo
            COMPLII(m) = (3.0d0 * PI * COMPI(m) * raio_cub) / (2.0d0 * MEY(m) * ESPHI(m))
            
            print '(A,I3)', "  Elemento ", m
            print '(A,ES14.6)', "    Resistencia:  ", RESISTI(m)
            print '(A,ES14.6)', "    Indutancia:   ", INDUCTI(m)
            print '(A,ES14.6)', "    Compliancia:  ", COMPLII(m)
            print *, ""
        end do

        print *, "========================================="
        print *, "Calculo das propriedades concluido!"
        print *, "========================================="
    end subroutine GEOMAT

    subroutine FORMAB(n, A, b, RESISTI, INDUCTI, COMPLII)
        implicit none
        integer :: n
        real(kind=8), dimension(n,n) :: A
        real(kind=8), dimension(n) :: b
        real(kind=8), dimension(:) :: RESISTI, INDUCTI, COMPLII
        integer :: i, j
        
        print *, ""
        print *, "========================================="
        print *, "MONTANDO MATRIZ A DO SISTEMA"
        print *, "========================================="
        
        ! Inicializa matriz A com zeros
        A = 0.0d0
        
        ! Diagonal principal = 1 (será sobrescrita pelos coeficientes específicos)
        do i = 1, n
            A(i,i) = 1.0d0
        end do

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

        
        ! Vetor b - inicialização
        b = 0.0d0
        if (n >= 1) b(1) = 1.0d0
        if (n >= 2) b(2) = 1.0d0

        print *, "Matriz A montada com sucesso!"
        print *, "Dimensao: ", n, "x", n
        print *, "========================================="
    end subroutine FORMAB

    subroutine MIVEULER(n, A, b, resultados_file)
        implicit none
        integer :: n
        real(kind=8) :: A(n,n)
        real(kind=8) :: b(n)
        character(len=*) :: resultados_file
        integer :: i, j, nt
        real(kind=8), allocatable :: Y(:), dYdt(:), Y_analitico(:), Y0(:)
        real(kind=8) :: tempo, tempo_total, dt
        
        ! Parâmetros
        tempo_total = 1.0d0
        dt = 0.1d0
        nt = int(tempo_total / dt)
        
        open(unit=20, file=resultados_file, status='replace', action='write')
        
        write(20, '(A)') '========================================='
        write(20, '(A)') 'RESULTADOS COMPLETOS DA SIMULACAO'
        write(20, '(A)') '========================================='
        write(20, '(A,I5)') 'Numero de elementos: ', n
        write(20, '(A,F10.6)') 'Tempo total da simulacao: ', tempo_total
        write(20, '(A,I10)') 'Numero de passos: ', nt
        write(20, '(A,F10.6)') 'Passo de tempo dt: ', dt
        write(20, '(A)') ''

        allocate(Y(n), dYdt(n), Y_analitico(n), Y0(n))

        print *, ""
        print *, "========================================="
        print *, "METODO DE EULER - PARAMETROS AUTOMATICOS"
        print *, "========================================="
        print *, "Tempo total: ", tempo_total, " segundos"
        print *, "Intervalo dt: ", dt, " segundos"
        print *, "Numero de passos: ", nt
        print *, ""
        
        Y0 = 0.0d0
        Y = Y0
        
        print *, "Vetor inicial Y0:"
        print '(5X, 100ES15.6)', (Y0(i), i = 1, n)

        write(20, '(A)') 'Condicao inicial (tempo = 0.0):'
        write(20, '(5X, 100ES15.6)') (Y0(i), i = 1, n)
        write(20, '(A)') ''
        
        write(20, '(A)') '========================================='
        write(20, '(A)') 'EVOLUCAO TEMPORAL - METODO DE EULER'
        write(20, '(A)') '========================================='
        write(20, '(A5, 2X, A10, 100(4X, A2, I3.3, A1))') &
              'Passo', 'Tempo', ('Y(', i, ')', i = 1, n)
        write(20, '(A)') '-----------------------------------------' 

        print *, ""
        print *, "Executando metodo de Euler..."
        print *, "----------------------------------------"
        
        tempo = 0.0d0
        do j = 0, nt
            write(20, '(I5, 2X, F10.6, 100(2X, ES15.6))') j, tempo, (Y(i), i = 1, n)
            print '(A,I5,A,F10.6)', "Passo ", j, " (t=", tempo, ")"
            
            if (j < nt) then
                dYdt = matmul(A, Y) + b
                Y = Y + dYdt * dt
                tempo = tempo + dt
            end if
        end do
        
        write(20, '(A)') ''
        write(20, '(A)') '========================================='
        write(20, '(A)') 'SOLUCAO ANALITICA'
        write(20, '(A)') '========================================='
        write(20, '(A)') 'Para o sistema dY/dt = Y + 1 com Y0 = 0:'
        write(20, '(A)') 'Solucao analitica: Y(x) = exp(x) - 1'
        write(20, '(A)') ''
        write(20, '(A10, 2X, A15, 2X, A15, 2X, A15)') 'Componente', 'Y0', 'Y_numerico', 'Y_analitico'
        write(20, '(A)') '-----------------------------------------------------------'
        
        print *, ""
        print *, "========================================="
        print *, "SOLUCAO ANALITICA NO TEMPO FINAL"
        print *, "========================================="
        
        do i = 1, n
            Y_analitico(i) = exp(tempo_total) - 1.0d0
            write(20, '(I10, 2X, ES15.6, 2X, ES15.6, 2X, ES15.6)') &
                  i, Y0(i), Y(i), Y_analitico(i)
        end do
        
        write(20, '(A)') ''
        write(20, '(A)') '========================================='
        write(20, '(A)') 'FIM DOS RESULTADOS'
        write(20, '(A)') '========================================='
        
        close(20)
        
        print *, ""
        print *, "Simulacao concluida!"
        print *, "Resultados salvos em '", trim(resultados_file), "'"
        print *, "========================================="

        deallocate(Y, dYdt, Y_analitico, Y0)
    end subroutine MIVEULER
 
    subroutine SALVAR_PROPRIEDADES(n, COMPI, RAIOI, ESPHI, MEY, RESISTI, INDUCTI, COMPLII, props_file)
        implicit none
        integer :: n
        real(kind=8), dimension(n) :: COMPI, RAIOI, ESPHI, MEY
        real(kind=8), dimension(n) :: RESISTI, INDUCTI, COMPLII
        character(len=*) :: props_file
        integer :: i
        
        open(unit=30, file=props_file, status='replace', action='write')
        
        write(30, '(A)') '========================================================================================='
        write(30, '(A)') 'PROPRIEDADES DOS ELEMENTOS'
        write(30, '(A)') '========================================================================================='
        write(30, '(A, I5)') 'Numero de elementos: ', n
        write(30, '(A)') ''
        write(30, '(A8, 7A12)') 'Elemento', 'COMP', 'RAIO', 'ESP', 'MEY', 'RESIST', 'INDUCT', 'COMPLII'
        write(30, '(A)') '-------------------------------------------------------------------------------------------'
        
        do i = 1, n
            write(30, '(I8, 7ES12.4)') i, COMPI(i), RAIOI(i), ESPHI(i), MEY(i), &
                                      RESISTI(i), INDUCTI(i), COMPLII(i)
        end do
        
        close(30)
        
        print *, ""
        print *, "Propriedades salvas em '", trim(props_file), "'"
    end subroutine SALVAR_PROPRIEDADES
 
    subroutine VISUALIZAR_MATRIZ(n, A, arquivo_viz)
        implicit none
        integer :: n
        real(kind=8), dimension(n,n) :: A
        character(len=*) :: arquivo_viz
        integer :: i, j, n_nao_zeros
        real(kind=8) :: tol
        
        tol = 1.0d-15
        n_nao_zeros = 0
        
        ! Contar elementos não-nulos
        do i = 1, n
            do j = 1, n
                if (abs(A(i,j)) > tol) then
                    n_nao_zeros = n_nao_zeros + 1
                endif
            end do
        end do
        
        open(unit=40, file=arquivo_viz, status='replace', action='write')
        
        write(40, '(A)') '================================================================================'
        write(40, '(A)') '                    VISUALIZACAO DA MATRIZ ESPARSA A'
        write(40, '(A)') '================================================================================'
        write(40, '(A)') ''
        write(40, '(A,I6,A,I6)') 'Dimensao da matriz: ', n, ' x ', n
        write(40, '(A,I8)') 'Total de elementos: ', n*n
        write(40, '(A,I8)') 'Elementos nao-nulos: ', n_nao_zeros
        write(40, '(A,F8.4,A)') 'Esparsidade: ', 100.0d0*(1.0d0 - real(n_nao_zeros,8)/real(n*n,8)), '%'
        write(40, '(A)') ''
        write(40, '(A)') '================================================================================'
        write(40, '(A)') 'LISTA DE ELEMENTOS NAO-NULOS (formato: LINHA, COLUNA, VALOR)'
        write(40, '(A)') '================================================================================'
        write(40, '(A)') ''
        write(40, '(A6, 2X, A6, 2X, A20)') 'Linha', 'Coluna', 'Valor'
        write(40, '(A)') '--------------------------------------------------------------------------------'
        
        do i = 1, n
            do j = 1, n
                if (abs(A(i,j)) > tol) then
                    write(40, '(I6, 2X, I6, 2X, ES20.10)') i, j, A(i,j)
                endif
            end do
        end do
        
        write(40, '(A)') ''
        write(40, '(A)') '================================================================================'
        write(40, '(A)') '                    MAPA DE ESPARSIDADE DA MATRIZ'
        write(40, '(A)') '================================================================================'
        write(40, '(A)') ''
        write(40, '(A)') 'Legenda: * = elemento nao-nulo,  . = zero'
        write(40, '(A)') ''
        
        ! Cabeçalho com números das colunas (dezenas)
        write(40, '(A6)', advance='no') '      '
        do j = 1, n
            if (mod(j, 10) == 0) then
                write(40, '(I1)', advance='no') mod(j/10, 10)
            else
                write(40, '(A1)', advance='no') ' '
            endif
        end do
        write(40, '(A)') ''
        
        ! Cabeçalho com números das colunas (unidades)
        write(40, '(A6)', advance='no') '      '
        do j = 1, n
            write(40, '(I1)', advance='no') mod(j, 10)
        end do
        write(40, '(A)') ''
        write(40, '(A)') ''
        
        ! Mapa completo da matriz (todas as linhas e colunas)
        do i = 1, n
            write(40, '(I5,A1)', advance='no') i, ':'
            do j = 1, n
                if (abs(A(i,j)) > tol) then
                    write(40, '(A1)', advance='no') '*'
                else
                    write(40, '(A1)', advance='no') '.'
                endif
            end do
            write(40, '(A)') ''
        end do
        
        write(40, '(A)') ''
        write(40, '(A)') '================================================================================'
        
        close(40)
        
        print *, ""
        print *, "Visualizacao da matriz salva em '", trim(arquivo_viz), "'"
        print *, "  Elementos nao-nulos: ", n_nao_zeros, " de ", n*n
        print *, "  Esparsidade: ", 100.0d0*(1.0d0 - real(n_nao_zeros,8)/real(n*n,8)), "%"
    end subroutine VISUALIZAR_MATRIZ
 
end module routines_p

program ZERODMOD
    use routines_p
    implicit none
    integer :: n, n_test
    real(kind=8) :: MI, RO
    real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
    real(kind=8), dimension(:), allocatable :: RESISTI, INDUCTI, COMPLII
    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8), dimension(:), allocatable :: b

    print *, ""
    print *, "========================================="
    print *, "  PROGRAMA ZERODMOD - INICIANDO"
    print *, "========================================="
    print *, ""

    ! 1. Ler dados no NOVO formato (6 valores por linha)
    call DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO)
    
    ! 2. Calcular propriedades
    allocate(RESISTI(n), INDUCTI(n), COMPLII(n))
    call GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, RESISTI, INDUCTI, COMPLII)
    
    ! 3. Salvar propriedades
    call SALVAR_PROPRIEDADES(n, COMPI, RAIOI, ESPHI, MEY, RESISTI, INDUCTI, COMPLII, 'propriedades_elementos.txt')
    
    ! 4. Definir dimensão para teste
    print *, ""
    print *, "========================================="
    print *, "MONTAGEM DO SISTEMA LINEAR"
    print *, "========================================="
    
    n_test = 122  ! Dimensão necessária para acomodar todos os índices
    
    ! 5. Alocar matriz A e vetor b
    allocate(A(n_test, n_test), b(n_test))
    
    ! 6. Montar sistema linear
    call FORMAB(n_test, A, b, RESISTI, INDUCTI, COMPLII)
    
    ! 7. Visualizar matriz esparsa
    call VISUALIZAR_MATRIZ(n_test, A, 'visualizacao_matriz.txt')
    
    ! 8. Executar método de Euler
    call MIVEULER(n_test, A, b, 'resultados_completos.txt')

    ! Liberar memória
    deallocate(COMPI, RAIOI, ESPHI, MEY)
    deallocate(RESISTI, INDUCTI, COMPLII)
    deallocate(A, b)

    print *, ""
    print *, "========================================="
    print *, "  EXECUCAO CONCLUIDA COM SUCESSO"
    print *, "========================================="
    print *, "Arquivos gerados:"
    print *, "  1. resultados_completos.txt"
    print *, "  2. propriedades_elementos.txt"
    print *, "  3. visualizacao_matriz.txt"
    print *, "========================================="
end program ZERODMOD
