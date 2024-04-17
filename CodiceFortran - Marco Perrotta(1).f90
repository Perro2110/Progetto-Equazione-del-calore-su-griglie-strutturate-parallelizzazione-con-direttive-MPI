program main
    ! ---------------------------------------------------------------------- !
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: i,j, n, timestep = 0, idummy ! utili
    INTEGER, PARAMETER :: IMAX=350    ! Numero di celle verticali/orizzontali || Estraibile dimMax
    INTEGER, PARAMETER :: NMAX=1e6    ! Numero max di time steps
    REAL    :: xL, xR  ! left and right coordinate del dominio
    REAL    :: yT, yB  ! top and bottom coordinate del dominio

    REAL    :: dx, dx2 ! mesh space (and its squared value)
    REAL    :: dy, dy2 ! mesh space (and its squared value)

    REAL    :: CFL     ! CFL number (<=1 for stability of EXPLICIT SCHEMES) 
    REAL    :: time    ! Tempo corrente
    REAL    :: dt      ! time step
    REAL    :: tend    ! Tempo finale
    REAL    :: TL, TR  ! "Condizioni di confine"
    REAL    :: kappa   ! coefficiente di conduzione del calore = k
    INTEGER :: NDIM = 2 ! quante dimensioni va diviso
    !
    REAL, ALLOCATABLE    :: T(:,:), Tnew(:,:) ! soluzioni correnti e nuove 
    REAL, ALLOCATABLE    :: x(:),y(:)         ! coordinate verticali e orizontali
    LOGICAL, ALLOCATABLE :: periods(:)
    INTEGER, ALLOCATABLE :: dims(:)
                                                           !! sezione MPI 
    INTEGER                   :: TCPU, BCPU, MsgLength, nMsg
    REAL, ALLOCATABLE         :: send_messageB(:), send_messageT(:), recv_messageB(:), recv_messageT(:)
    INTEGER                   :: send_request(2), recv_request(2)  ! tag dei messaggi che sono stati mandati
    INTEGER                   :: send_status_list(MPI_STATUS_SIZE,2), recv_status_list(MPI_STATUS_SIZE,2) 
    INTEGER                   :: source
    INTEGER,DIMENSION(8)      :: timeI,timeF
    


    TYPE tMPI                           ! STRUTTURA
        INTEGER :: myrank, nCPU, iErr
        INTEGER :: iMMax
        INTEGER :: AUTO_REAL
        INTEGER :: nElem, iStart, iEnd  ! numero di celle che ha ogni cpu,
        !
        INTEGER, ALLOCATABLE :: myCoords (:) ! coordinate cella in cui si trova il temp
    END TYPE tMPI

    TYPE(tMPI) :: MPI
    REAL       :: realtest
    INTEGER    :: COMM_CART !! CANALE DI COMUNICAZIONE PER IL PIANO CARTESIANO

    ! ---------------------------------------------------------------------- !
    
    !!
    ! CLASSICAL MPI INITIALIZATION
    !!
    
    CALL MPI_INIT(MPI%iErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPI%myrank,MPI%iErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPI%nCPU,MPI%iErr)

    ! per capire il "TIPO" DI messaggi da inviare
    ! per rendere generale la compilazione del codice
    !_______________________________________________________
    SELECT CASE(KIND(realtest))
    CASE(4)
        MPI%AUTO_REAL = MPI_REAL
    CASE(8)
        MPI%AUTO_REAL = MPI_DOUBLE_PRECISION
    END SELECT
    ! ______________________________________________________

    CALL DATE_AND_TIME(VALUES=timeI)
    !!
    ! DOMAIN DECOMPOSITION
    !!
    
    IF(MOD(IMAX,MPI%nCPU).NE.0) THEN  ! vogliamo che il numero di celle sia multiplo del numero di cpu 
        PRINT *, ' ERROR. Number of cells (IMAX) must be a multiple of number of CPU (nCPU)!'
        CALL MPI_FINALIZE(MPI%iErr)
        STOP
    ELSE
        IF(MPI%myrank.EQ.0) WRITE(*,*) ' Parallel simulation with MPI. '
    ENDIF

    CONTINUE
    
    !!ALLOCO
    ALLOCATE(dims(NDIM),periods(NDIM),MPI%myCoords(NDIM))
    !!ASSEGNO VALORI
    dims    = (/1,MPI%nCpu/)
    periods = (/.FALSE.,.TRUE./)
    
    continue
    !!
    CALL MPI_CART_CREATE(MPI_COMM_WORLD,NDIM,dims,periods,.TRUE.,COMM_CART,MPI%iErr)
    CALL MPI_COMM_RANK(COMM_CART,MPI%myrank,MPI%iErr)
    !

    ! TROVA I VICINI (DINGLEBERG)
    CALL MPI_CART_SHIFT(COMM_CART,1,1,source,TCPU,MPI%iErr)
    CALL MPI_CART_SHIFT(COMM_CART,1,-1,source,BCPU,MPI%iErr)
    !

    !FARE COORD
    CALL MPI_CART_COORDS(COMM_CART,MPI%myrank,NDIM,MPI%myCoords,MPI%iErr)
    MPI%nElem  = IMAX                           ! orizontalmente va dall' inizio alla fine
    MPI%iMMax  = IMAX/MPI%nCpu                  ! verticalmente  numero di celle diviso per ogni cpu
   
    MPI%iStart = 1 + MPI%myCoords(1)*MPI%iMMax  ! sempre verticalmente
    MPI%iEnd   = MPI%iStart + MPI%iMMax - 1     !

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MsgLength = IMAX !! NEL CASO MPI%NUMELEM
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CONTINUE
    CFL = 0.9 ! CFL number

    ! Domain definition
    xL   = 0.
    xR   = 2.
    
    yT = 0.
    yB = 2.
    ! Boundary conditions
    TL   = 100.
    TR   = 50.

    time = 0.0
    tend = 0.05

    kappa= 1.0

    ALLOCATE( x(IMAX) , y(IMAX) )                       ! GRIGLIA
    ALLOCATE( Tnew(IMAX,IMAX) )                         ! SOL
    ALLOCATE( T(MPI%iStart-1:MPI%iEnd+1, 1:MPI%nElem) ) ! SOL
    ALLOCATE(send_messageT(MPI%nElem), send_messageB(MPI%nElem))
    ALLOCATE(recv_messageT(MPI%nElem), recv_messageB(MPI%nElem))

    IF(MPI%myrank.EQ.0) WRITE(*,*) ' Building the computational domain... ' !debug!
    CONTINUE
    ! FAMO I DELTA
    dx = (xR-xL)/REAL(IMAX-1)
    dx2= dx**2

    dy = (yB-yT)/REAL(IMAX-1) 
    dy2= dy**2

    x(1) = xL
    DO i = 1, IMAX-1
        x(i+1) = x(i) + dx
    ENDDO

    y(1) = yT
    DO i = 1, IMAX-1
        y(i+1) = y(i) + dy
    ENDDO

    IF(MPI%myrank.EQ.0) WRITE(*,*) ' Assigning the initial condition... '

    DO j = 1, MPI%nElem
        IF(x(j).LE.1) THEN
            T(:,j) = TL
        ELSE
            T(:,j) = TR
        ENDIF
    ENDDO

    DO j = 1, IMAX
        CALL MPI_ALLGATHER( T(:,j), MPI%iMMax, MPI%AUTO_REAL, Tnew(:,j), MPI%iMMax, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr )
    ENDDO
    IF(MPI%myrank.EQ.0) CALL PlotOutputMPI(IMAX,x,y,Tnew,timestep,MPI%myrank,1,IMAX,MPI%nCpu)                                                         !!
    

    IF(MPI%myrank.EQ.0) WRITE(*,*) ' START OF THE COMPUTATION '

    CONTINUE
    DO n = 1, NMAX   ! main loop in time

        IF(time.GE.tend) EXIT

        dt = 0.5*CFL*dx2/kappa 
        IF(time+dt.GT.tend) THEN
            dt = tend-time  ! adjust the last time step in order to exactly match tend
        ENDIF

        ! NUMERICAL SCHEME: EXPLICIT finite difference method

        !!
        ! Message exchange
        !! 
        
        ! 1) SETUP ALL COMMUNICATIONS
        ! CREIAMO PRIMA I MESSAGGI
        send_messageT = T(MPI%iStart, 1 : MPI%nElem)
        send_messageB = T(MPI%iEnd,   1 : MPI%nElem)

        ! no blocking comunication
        CALL MPI_ISEND(send_messageT,MsgLength,MPI%AUTO_REAL,TCPU,1,COMM_CART,send_request(1),MPI%iErr) ! BOTTER AL TOPPER
        CALL MPI_ISEND(send_messageB,MsgLength,MPI%AUTO_REAL,BCPU,2,COMM_CART,send_request(2),MPI%iErr) ! viceversa

        CALL MPI_IRECV(recv_messageT,MsgLength,MPI%AUTO_REAL,TCPU,2,COMM_CART,recv_request(1),MPI%iErr)
        CALL MPI_IRECV(recv_messageB,MsgLength,MPI%AUTO_REAL,BCPU,1,COMM_CART,recv_request(2),MPI%iErr)

        ! 2) USEFUL TIME! In the meantime, we update all internal (non MPI) elements
        CONTINUE
        DO i = MPI%iStart+1, MPI%iEnd-1
            DO j = 2, MPI%nElem-1
                Tnew(i,j) = T(i,j) + kappa*dt/dx2 * (T(i+1,j) - 2.*T(i,j) + T(i-1,j)) + kappa*dt/dy2 * (T(i,j+1) - 2.*T(i,j) + T(i,j-1))
            ENDDO
        ENDDO

        ! 3) UPDATE MPI boundary elements
        !    aspettiamo che la comunicazione sia finita
        
        nMsg = 2
        CALL MPI_WAITALL(nMsg,send_request,send_status_list,MPI%iErr)
        CALL MPI_WAITALL(nMsg,recv_request,recv_status_list,MPI%iErr)

        ! AGGIORNAMENTO BORDI

        T(MPI%iStart-1,:) = recv_messageT
        T(MPI%iEnd+1,:)   = recv_messageB

        ! MANCANO SOLO I BORDI "i = 1 e i = 5"
            ! fisici
        Tnew(:,1)         = TL
        Tnew(:,MPI%nElem) = TR
            ! non
        DO i = MPI%iStart, MPI%iEnd,(MPI%iEnd - MPI%iStart)
            DO j = 2, MPI%nElem-1
                Tnew(i,j) = T(i,j) + kappa*dt/dx2 * (T(i+1,j) - 2.*T(i,j) + T(i-1,j)) + kappa*dt/dy2 * (T(i,j+1) - 2.*T(i,j) + T(i,j-1))
            ENDDO
        ENDDO

        ! Update time and current solution
        time = time + dt
        T(MPI%istart:MPI%iend,:) = Tnew(MPI%istart:MPI%iend,:) ! sovrascrivo la nuova soluzione
        timestep = n
        PRINT *,n," tempo: ",time !! per vedere se va avanti
    ENDDO !n

    IF(MPI%myrank.EQ.0) WRITE(*,*) ' END OF THE COMPUTATION '
    
    ! colleziono i dati dai vari domini
    DO j = 1, IMAX
        CALL MPI_ALLGATHER( T(:,j), MPI%iMMax, MPI%AUTO_REAL, Tnew(:,j), MPI%iMMax, MPI%AUTO_REAL, MPI_COMM_WORLD, MPI%iErr )
    ENDDO
    CALL DATE_AND_TIME(VALUES=timeF)
    
    timestep = ( timeF(6)*60000 + timeF(7)*1000 + timeF(8) ) - ( timeI(6)*60000 + timeI(7)*1000 + timeI(8) ) 
    
    IF(MPI%myrank.EQ.0) CALL PlotOutputMPI(IMAX,x,y,Tnew,timestep,MPI%myrank,1,IMAX,MPI%nCpu)
    

    ! Empty memory
    DEALLOCATE( T, Tnew, x , y )

    CALL MPI_FINALIZE(MPI%iErr)
end program main
    

SUBROUTINE PlotOutputMPI(IMAX,x,y,T,n,myrank,istart,iend,ncpu)
    IMPLICIT NONE
    INTEGER             :: IMAX, i, j, n, myrank, istart, iend, OutUnit,ncpu
    REAL                :: x(IMAX),y(IMAX), T(istart:iend,IMAX)
    CHARACTER(LEN=10)   :: citer, cmyrank
    CHARACTER(LEN=200)  :: IOFileName

    WRITE(citer,'(I4.4)') n
    WRITE(cmyrank,'(I4.4)') myrank
    IOFileName = 'Heat2D_output-'//TRIM(citer)//'-'//TRIM(cmyrank)//'.dat'

    OutUnit = 100+myrank
    OPEN(UNIT=OutUnit, FILE=IOFileName, STATUS='unknown', ACTION='write')

    WRITE(OutUnit,*) iend-istart+1
    WRITE(OutUnit,*) n
    WRITE(OutUnit,*) ncpu
    
    DO i = 1, IMAX
        WRITE(OutUnit,*) x(i)
    ENDDO
    
    DO i = 1, IMAX
        WRITE(OutUnit,*) y(i)
    ENDDO
    
    DO i = 1, IMAX
        DO j = 1, IMAX
            WRITE(OutUnit,*) T(i,j)
        ENDDO
    ENDDO

END SUBROUTINE PlotOutputMPI