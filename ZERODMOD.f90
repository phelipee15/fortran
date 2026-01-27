module routines_p
    implicit none
    public :: DADOSMOD, GEOMAT, FORMAB, MIVEULER, SALVAR_PROPRIEDADES
    private
    real(kind=8), parameter :: PI = 3.141592653589793d0
    
contains

    subroutine DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO
        real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
        integer :: i, ios
        character(len=100) :: linha

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

        ! Lê dados geométricos (comprimento, raio, espessura)
        print *, "Lendo dados geometricos..."
        do i = 1, n
            read(10,*,iostat=ios) COMPI(i), RAIOI(i), ESPHI(i)
            if (ios /= 0) then
                print *, "ERRO: Falha ao ler linha ", i, " dos dados geometricos"
                stop
            endif
            print '(A,I3,A,F8.3,A,F8.3,A,F8.3)', "  Elemento ", i, &
                  ": COMP=", COMPI(i), " RAIO=", RAIOI(i), " ESPH=", ESPHI(i)
        end do

        ! Lê propriedades do material (Módulo de Young para cada elemento)
        print *, "Lendo modulo de Young para cada elemento..."
        do i = 1, n
            read(10,*,iostat=ios) MEY(i)
            if (ios /= 0) then
                print *, "ERRO: Nao foi possivel ler MEY para o elemento ", i
                stop
            endif
            print '(A,I3,A,F12.3)', "  Elemento ", i, ": MEY = ", MEY(i)
        end do

        ! Lê Permeabilidade magnética e Densidade (constantes)
        read(10,*,iostat=ios) MI, RO
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel ler MI e RO"
            stop
        endif
        
        close(10)

        print *, "========================================="
        print *, "Propriedades do Material:"
        print *, "  MI:", MI
        print *, "  RO:", RO
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
            
            ! Resistência
            RESISTI(m) = (8.0d0 * MI * COMPI(m)) / (PI * raio_quarta)
            
            ! Indutância
            INDUCTI(m) = (9.0d0 * RO * COMPI(m)) / (4.0d0 * PI * raio_quad)
            
            ! Compliância (flexibilidade)
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
        real(kind=8), dimension(n) :: RESISTI, INDUCTI, COMPLII
        integer :: i, j
        
        ! Construção do problema hipotético: A = Identidade, b pré-definido
        A = 0.0d0
        do i = 1, n
            A(i,i) = 1.0d0  ! Diagonal principal = 1

        ! COEFICIENTES NÃO NULOS DA MATRIZ A
        A(1,4) = 1/COMPLII(1)
        A(1,10) = 1/COMPLII(1)
        A(1,16) = 1/COMPLII(1)
        A(2,1) = 1/INDUCTI(1)
        A(2,4) = -RESISTI(1)/INDUCTI(1)
        A(2,10) = -RESISTI(1)/INDUCTI(1)
        A(2,16) = -RESISTI(1)/INDUCTI(1)
        A(3,4) = 1/COMPLII(2)
        A(3,6) = -1/COMPLII(2)
        A(3,8) = -1/COMPLII(2)
        A(4,1) = 1/INDUCTI(2)
        A(4,3) = -1/INDUCTI(2)
        A(4,4) = -RESISTI(2)/INDUCTI(2)
        A(5,6) = 1/COMPLII(3)
        A(6,5) = 1/INDUCTI(3)
        A(6,6) = -RESISTI(3)/INDUCTI(3)
        A(7,8) = 1/COMPLII(4)
        A(8,7) = 1/INDUCTI(4)
        A(8,8) = -RESISTI(4)/INDUCTI(4)
        A(9,10) = 1/COMPLII(5)
        A(9,12) = -1/COMPLII(5)
        A(9,14) = -1/COMPLII(5)
        A(10,1) = -1/INDUCTI(5)
        A(10,9) = -1/INDUCTI(5)
        A(10,10) = -RESISTI(5)/INDUCTI(5)
        A(11,12) = 1/COMPLII(6)
        A(12,11) = 1/INDUCTI(6)
        A(12,12) = -RESISTI(6)/INDUCTI(6)
        A(13,14) = 1/COMPLII(7)
        A(14,13) = 1/INDUCTI(7)
        A(14,14) = -RESISTI(7)/INDUCTI(7)
        A(15,16) = 1/COMPLII(8)
        A(15,18) = -1/COMPLII(8)
        A(15,20) = -1/COMPLII(8)
        A(16,1) = 1/INDUCTI(8)
        A(16,15) = -1/INDUCTI(8)
        A(16,18) = -RESISTI(8)/INDUCTI(8)
        A(16,20) = -RESISTI(8)/INDUCTI(8)
        A(17,18) = 1/COMPLII(9)
        A(17,42) = -1/COMPLII(9)
        A(17,44) = -1/COMPLII(9)
        A(18,15) = 1/INDUCTI(9)
        A(18,17) = -1/INDUCTI(9)
        A(18,42) = -RESISTI(9)/INDUCTI(9)
        A(18,44) = -RESISTI(9)/INDUCTI(9)
        A(19,20) = 1/COMPLII(10)
        A(19,22) = -1/COMPLII(10)
        A(19,24) = -1/COMPLII(10)
        A(20,15) = 1/INDUCTI(10)
        A(20,19) = -1/INDUCTI(10)
        A(20,22) = -RESISTI(10)/INDUCTI(10)
        A(20,24) = -RESISTI(10)/INDUCTI(10)
        A(21,22) = 1/COMPLII(11)
        A(21,26) = -1/COMPLII(11)
        A(21,28) = -1/COMPLII(11)
        A(22,19) = 1/INDUCTI(11)
        A(22,21) = -1/INDUCTI(11)
        A(22,26) = -RESISTI(11)/INDUCTI(11)
        A(22,28) = -RESISTI(11)/INDUCTI(11)
        A(23,24) = 1/COMPLII(12)
        A(23,40) = -1/COMPLII(12)
        A(23,38) = -1/COMPLII(12)
        A(24,19) = 1/INDUCTI(12)
        A(24,23) = -1/INDUCTI(12)
        A(24,38) = -RESISTI(12)/INDUCTI(12)
        A(24,40) = -RESISTI(12)/INDUCTI(12)
        A(25,26) = -1/COMPLII(13)
        A(26,25) = 1/INDUCTI(13)
        A(26,26) = -RESISTI(13)/INDUCTI(13)
        A(27,28) = 1/COMPLII(14)
        A(27,30) = -1/COMPLII(14)
        A(27,32) = -1/COMPLII(14)
        A(28,21) = 1/INDUCTI(14)
        A(28,27) = -1/INDUCTI(14)
        A(28,30) = -RESISTI(14)/INDUCTI(14)
        A(28,32) = -RESISTI(14)/INDUCTI(14)
        A(29,30) = 1/COMPLI(15)
        A(30,29) = 1/INDUCTI(15)
        A(30,30) = -RESISTI(15)/INDUCTI(15)
        A(31,32) = 1/COMPLI(16)
        A(31,34) = -1/COMPLI(16)
        A(31,36) = -1/COMPLI(16)
        A(32,27) = 1/INDUCTI(16)
        A(32,31) = -1/INDUCTI(16)
        A(32,34) = -RESISTI(16)/INDUCTI(16)
        A(32,36) = -RESISTI(16)/INDUCTI(16)
        A(33,34) = 1/COMPLI(17)
        A(34,33) = 1/INDUCTI(17)
        A(34,34) = -RESISTI(17)/INDUCTI(17)
        A(35,36) = 1/COMPLI(18)
        A(36,35) = 1/INDUCTI(18)
        A(36,36) = -RESISTI(18)/INDUCTI(18)
        A(37,38) = 1/COMPLI(19)
        A(38,37) = 1/INDUCTI(19)
        A(38,38) = -RESISTI(19)/INDUCTI(19)
        A(39,40) = 1/COMPLI(20)
        A(40,39) = 1/INDUCTI(20)
        A(40,40) = -RESISTI(20)/INDUCTI(20)
        A(41,42) = 1/COMPLI(21)
        A(41,50) = -1/COMPLI(21)
        A(41,52) = -1/COMPLI(21)
        A(42,17) = 1/INDUCTI(21)
        A(42,41) = -1/INDUCTI(21)
        A(42,50) = -RESISTI(21)/INDUCTI(21)
        A(42,52) = -RESISTI(21)/INDUCTI(21)
        A(43,44) = 1/COMPLI(22)
        A(43,46) = -1/COMPLI(22)
        A(43,48) = -1/COMPLI(22)
        A(44,17) = 1/INDUCTI(22)
        A(44,43) = -1/INDUCTI(22)
        A(44,46) = -RESISTI(22)/INDUCTI(22)
        A(44,48) = -RESISTI(22)/INDUCTI(22)
        A(45,46) = 1/COMPLI(23)
        A(46,45) = 1/INDUCTI(23)
        A(46,46) = -RESISTI(23/INDUCTI(23)
        A(47,48) = 1/COMPLI(24)
        A(48,47) = 1/INDUCTI(24)
        A(48,48) = -RESISTI(24)/INDUCTI(24)
        A(49,50) = 1/COMPLI(25)
        A(49,68) = -1/COMPLI(25)
        A(49,66) = -1/COMPLI(25)
        A(50,41) = 1/INDUCTI(25)
        A(50,49) = -1/INDUCTI(25
        A(50,66) = -RESISTI(25)/INDUCTI(25)
        A(50,68) = -RESISTI(25)/INDUCTI(25)
        A(51,52) = 1/COMPLI(26)
        A(51,54) = -1/COMPLI(26)
        A(51,56) = -1/COMPLI(26)
        A(52,48) = 1/INDUCTI(26)
        A(52,51) = -1/INDUCTI(26)
        A(52,54) = -RESISTI(26)/INDUCTI(26)
        A(52,56)  = -RESISTI(26)/INDUCTI(26)


        ! Vetor b pré-definido (valores específicos para cada elemento)
        if (n >= 1) b(1) = 1.0d0
        if (n >= 2) b(2) = 1.0d0
        ! Adicione mais valores se necessário para n > 2

        print *, ""
        print *, "========================================="
        print *, "MATRIZ A E VETOR b (PROBLEMA HIPOTÉTICO)"
        print *, "========================================="
        print *, "Matriz A (Identidade ", n, "x", n, "):"
        do i = 1, n
            print '(*(F8.3,1X))', (A(i,j), j = 1, n)
        end do

        print *, ""
        print *, "Vetor b (pré-definido):"
        print '(*(F8.3,1X))', (b(i), i = 1, n)
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
        
        ! Parâmetros pré-definidos
        tempo_total = 1.0d0  ! 1 segundo
        dt = 0.1d0           ! Intervalo de 0.1 segundo
        nt = int(tempo_total / dt)
        
        ! Abrir arquivo único para salvar TODOS os resultados
        open(unit=20, file=resultados_file, status='replace', action='write')
        
        ! ============================================================
        ! SEÇÃO 1: CABEÇALHO E INFORMAÇÕES DA SIMULAÇÃO
        ! ============================================================
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
        
        ! Vetor inicial Y0 pré-definido como zeros
        Y0 = 0.0d0
        Y = Y0
        
        print *, "Vetor inicial Y0 (pre-definido como zeros):"
        print '(5X, 100ES15.6)', (Y0(i), i = 1, n)

        ! Escrever condição inicial
        write(20, '(A)') 'Condicao inicial (tempo = 0.0) [pre-definido]:'
        write(20, '(5X, 100ES15.6)') (Y0(i), i = 1, n)
        write(20, '(A)') ''
        
        ! ============================================================
        ! SEÇÃO 2: EVOLUÇÃO TEMPORAL (MÉTODO DE EULER)
        ! ============================================================
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
        do j = 0, nt  ! Inclui passo 0 (condição inicial)
            ! Escrever no arquivo
            write(20, '(I5, 2X, F10.6, 100(2X, ES15.6))') j, tempo, (Y(i), i = 1, n)
            
            ! Imprimir na tela
            print '(A,I5,A,F10.6,A,100(ES13.5,1X))', "Passo ", j, " (t=", tempo, "): Y = ", Y
            
            ! Calcular próximo passo (exceto no último)
            if (j < nt) then
                dYdt = matmul(A, Y) + b
                Y = Y + dYdt * dt
                tempo = tempo + dt
            end if
        end do
        
        write(20, '(A)') ''
        
        ! ============================================================
        ! SEÇÃO 3: SOLUÇÃO ANALÍTICA
        ! ============================================================
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
        print *, "Y(x) = exp(x) - 1"
        print *, ""
        
        ! Calcular solução analítica para cada componente
        do i = 1, n
            Y_analitico(i) = exp(tempo_total) - 1.0d0
            write(20, '(I10, 2X, ES15.6, 2X, ES15.6, 2X, ES15.6)') &
                  i, Y0(i), Y(i), Y_analitico(i)
            print '(A,I3,A,ES15.6,A,ES15.6,A,ES15.6)', "Comp ", i, &
                  ": Y0=", Y0(i), " Y_num=", Y(i), " Y_ana=", Y_analitico(i)
        end do
        
        write(20, '(A)') ''
        write(20, '(A)') '========================================='
        write(20, '(A)') 'FIM DOS RESULTADOS'
        write(20, '(A)') '========================================='
        
        close(20)
        
        print *, ""
        print *, "----------------------------------------"
        print *, "Simulacao concluida!"
        print *, "Todos os resultados salvos em '", trim(resultados_file), "'"
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
        write(30, '(A8, 7A12)') 'Elemento', 'COMP', 'RAIO', 'ESP', 'MEY', 'RESIST', 'INDUCT', 'COMPLI'
        write(30, '(A)') '-------------------------------------------------------------------------------------------'
        
        do i = 1, n
            write(30, '(I8, 7ES12.4)') i, COMPI(i), RAIOI(i), ESPHI(i), MEY(i), &
                                      RESISTI(i), INDUCTI(i), COMPLII(i)
        end do
        
        close(30)
        
        print *, ""
        print *, "Propriedades salvas em '", trim(props_file), "'"
    end subroutine SALVAR_PROPRIEDADES
 
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

    ! 1. Ler dados do arquivo
    call DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO)
    
    ! 2. Calcular propriedades dos elementos
    allocate(RESISTI(n), INDUCTI(n), COMPLII(n))
    call GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, RESISTI, INDUCTI, COMPLII)
    
    ! 3. Salvar propriedades em arquivo
    call SALVAR_PROPRIEDADES(n, COMPI, RAIOI, ESPHI, MEY, RESISTI, INDUCTI, COMPLII, 'propriedades_elementos.txt')
    
    ! 4. FORÇAR n=2 para teste do problema hipotético
    print *, ""
    print *, "========================================="
    print *, "TESTE COM MATRIZ 2x2 (FORÇADO)"
    print *, "========================================="
    
    ! Definir n=2 para o teste
    n_test = 4
    
    ! 5. Alocar matriz A e vetor b para o teste
    allocate(A(n_test, n_test), b(n_test))
    
    ! 6. Montar sistema linear 2x2 (problema hipotético)
    call FORMAB(n_test, A, b, RESISTI, INDUCTI, COMPLII)
    
    ! 7. Executar método de Euler
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
    print *, "  1. resultados_completos.txt (TODOS OS RESULTADOS)"
    print *, "  2. propriedades_elementos.txt"
end program ZERODMOD



