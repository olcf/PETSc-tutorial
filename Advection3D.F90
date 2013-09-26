PROGRAM Advection3D

  USE RealKindModule, ONLY: &
    Zero, &
    One, &
    Two, &
    Three, &
    Pi, &
    Tiny, &
    Huge
  USE UniverseModule, ONLY: &
    UniverseType, &
    CreateUniverse, &
    DestroyUniverse
  USE MeshModule, ONLY: &
    MeshType, &
    X1, X2, X3
  USE FieldsModule, ONLY: &
    FieldsType, &
    D, V1, V2, V3
  USE AdvectionModule, ONLY: &
    AdvectExchangeCells3D, &
    AdvectNonExchangeCells3D

  IMPLICIT NONE

#include "finclude/petsc.h90"

  PetscBool :: &
    Done
  PetscErrorCode :: &
    Error    
  PetscInt :: &
    iX1, iX2, iX3, &
    iCycle, &
    MaxCycles, &
    X1B, X1E, &
    X2B, X2E, &
    X3B, X3E, &
    iField
  PetscInt, DIMENSION(3) :: &
    nZones
  PetscReal :: &
    MyTimeStep, &
    TimeStep, &
    Time, &
    EndTime, &
    MyL1Norm, &
    L1Norm
  PetscReal, DIMENSION(:,:,:,:), POINTER :: &
    InnerEdge, &
    Center, &
    OuterEdge, &
    FieldsOld, &
    FieldsNew
  Vec :: &
    InitialCondition, &
    FieldsScratch
  TYPE(MeshType), POINTER :: &
    M
  TYPE(FieldsType), POINTER :: &
    F
  TYPE(UniverseType), POINTER :: &
    U => NULL()

  !  Program solves dD/dt+Vx*dD/dx+Vy*dD/dy+Vz*dD/dz=0 on periodic unit domain 
  !    with constant advection velocity (V1, V2, V3) = (1, 0, 0), 
  !    and initial condition D = Sin(2 * Pi * X1)

  !  Create Mesh and Data Structures to Hold Physical Fields:

  nZones(1:3) = 256

  CALL CreateUniverse( &
         U, Name = 'Advection3D',          &
         nZones                            &
           = nZones,                       &
         InnerBoundaries                   &
           = (/ Zero, Zero, Zero /),       &
         OuterBoundaries                   &
           = (/ One,  One,  One  /),       &
         BoundaryConditions                &
           = (/ DMDA_BOUNDARY_PERIODIC,    &
                DMDA_BOUNDARY_PERIODIC,    &
                DMDA_BOUNDARY_PERIODIC /))

  M => U % Mesh   ! MeshType   Defined in MeshModule.F90
  F => U % Fields ! FieldsType Defined in FieldsModule.F90

  !  Create Shorter Names with Associate:

  ASSOCIATE( iBX1 => M % Positions % iBX(1), &
             iBX2 => M % Positions % iBX(2), &
             iBX3 => M % Positions % iBX(3), &
             iEX1 => M % Positions % iEX(1), &
             iEX2 => M % Positions % iEX(2), &
             iEX3 => M % Positions % iEX(3), &
             SW   => M % StencilWidth,       &
             DOF  => F % DegreesOfFreedom )

  !  Create Local Scratch Vector:

  CALL VecDuplicate(F % FieldsLocal, FieldsScratch, Error)

  !  Get Access to Mesh Coordinates:

  CALL DMDAVecGetArrayF90( &
         M % CoordinateDA, M % Positions % InnerEdgeGlobal, InnerEdge, Error)
  CALL DMDAVecGetArrayF90( &
         M % CoordinateDA, M % Positions % CenterGlobal,    Center,    Error)
  CALL DMDAVecGetArrayF90( &
         M % CoordinateDA, M % Positions % OuterEdgeGlobal, OuterEdge, Error)

  !  Get Access to Non Ghosted Fields Data:

  CALL DMDAVecGetArrayF90( &
         F % FieldsDA, F % FieldsGlobal, FieldsOld, Error)

  !  Set Initial Condition and Compute Time Step: 

  MyTimeStep = Huge

  DO iX3 = iBX3, iEX3
    DO iX2 = iBX2, iEX2
      DO iX1 = iBX1, iEX1

        FieldsOld(D,  iX1, iX2, iX3) &
          = SIN( Two * Pi * Center(X1, iX1, iX2, iX3) )

        FieldsOld(V1, iX1, iX2, iX3) &
          = One
        FieldsOld(V2, iX1, iX2, iX3) &
          = Zero
        FieldsOld(V3, iX1, iX2, iX3) &
          = Zero

        MyTimeStep &
          = MIN(MyTimeStep, &
                ( OuterEdge(X1, iX1, iX2, iX3) &
                  - InnerEdge(X1, iX1, iX2, iX3) ) &
                / ( Tiny + FieldsOld(V1, iX1, iX2, iX3) ) )

        MyTimeStep &
          = MIN(MyTimeStep, &
                ( OuterEdge(X2, iX1, iX2, iX3) &
                  - InnerEdge(X2, iX1, iX2, iX3) ) &
                / ( Tiny + FieldsOld(V2, iX1, iX2, iX3) ) )

        MyTimeStep &
          = MIN(MyTimeStep, &
                ( OuterEdge(X3, iX1, iX2, iX3) &
                  - InnerEdge(X3, iX1, iX2, iX3) ) &
                / ( Tiny + FieldsOld(V3, iX1, iX2, iX3) ) )

      END DO
    END DO
  END DO

  !  Release Access to Non-Ghosted Fields:

  CALL DMDAVecRestoreArrayF90( &
         F % FieldsDA, F % FieldsGlobal, FieldsOld, Error)

  !  Reduce Global Time Step:

  CALL MPI_ALLREDUCE( &
         MyTimeStep, TimeStep, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
         U % Communicator, Error)

  !  Populate Ghosted Fields with Initial Condition:

  CALL DMGlobalToLocalBegin( &
         F % FieldsDA, F % FieldsGlobal, INSERT_VALUES, F % FieldsLocal, Error)
  CALL DMGlobalToLocalEnd( &
         F % FieldsDA, F % FieldsGlobal, INSERT_VALUES, F % FieldsLocal, Error)

  !  Keep Copy of Initial Condition:

  CALL VecDuplicate(F % FieldsLocal, InitialCondition, Error)
  CALL VecCopy(     F % FieldsLocal, InitialCondition, Error)

  !  Begin Time Step Loop:

  Done      = .FALSE.
  iCycle    = 0
  MaxCycles = 100000
  Time      = Zero
  EndTime   = One

  DO WHILE (.NOT.Done)

    iCycle = iCycle + 1

    IF(U % MPI_Rank == 0)THEN
      PRINT*
      PRINT*, "    Cycle     = ", iCycle
      PRINT*, "    Time Step = ", TimeStep
      PRINT*, "    Time      = ", Time
    END IF

    !  Copy Ghosted Data to Local Work Vector:

    CALL VecCopy(F % FieldsLocal, FieldsScratch, Error)

    !  Get Access to Fields Data:

    CALL DMDAVecGetArrayF90( &
           F % FieldsDA, FieldsScratch,    FieldsOld, Error)
    CALL DMDAVecGetArrayF90( &
           F % FieldsDA, F % FieldsGlobal, FieldsNew, Error)

    !  Update Cells Whose Data is Communicated:

    CALL AdvectExchangeCells3D(FieldsNew, FieldsOld, InnerEdge, OuterEdge, &
           TimeStep, iBX1, iBX2, iBX3, iEX1, iEX2, iEX3, SW)

    !  Begin Ghost Exchange:

    CALL DMGlobalToLocalBegin( &
           F % FieldsDA, F % FieldsGlobal, INSERT_VALUES, F % FieldsLocal, &
           Error)

    !  Update Cells Whose Data is Not Comminucated:

    CALL AdvectNonExchangeCells3D(FieldsNew, FieldsOld, InnerEdge, OuterEdge, &
           TimeStep, iBX1, iBX2, iBX3, iEX1, iEX2, iEX3, SW)

    !  Release Access to Fields Data:

    CALL DMDAVecRestoreArrayF90( &
           F % FieldsDA, FieldsScratch,    FieldsOld, Error)
    CALL DMDAVecRestoreArrayF90( &
           F % FieldsDA, F % FieldsGlobal, FieldsNew, Error)

    !  End Ghost Exchange:

    CALL DMGlobalToLocalEnd( &
           F % FieldsDA, F % FieldsGlobal, INSERT_VALUES, F % FieldsLocal, &
           Error)

    !  Copy New to Old:

    CALL DMDAVecGetArrayF90( &
           F % FieldsDA, F % FieldsLocal,  FieldsOld, Error)
    CALL DMDAVecGetArrayF90( &
           F % FieldsDA, F % FieldsGlobal, FieldsNew, Error)

    DO iX3 = iBX3, iEX3
      DO iX2 = iBX2, iEX2
        DO iX1 = iBX1, iEX1
          DO iField = 0, DOF - 1

            FieldsOld(iField, iX1, iX2, iX3) &
              = FieldsNew(iField, iX1, iX2, iX3)

          END DO
        END DO
      END DO
    END DO

    CALL DMDAVecRestoreArrayF90( &
           F % FieldsDA, F % FieldsLocal,  FieldsOld, Error)
    CALL DMDAVecRestoreArrayF90( &
           F % FieldsDA, F % FieldsGlobal, FieldsNew, Error)

    Time = Time + TimeStep

    IF(iCycle >= MaxCycles .OR. Time >= EndTime) Done = .TRUE.

  END DO

  !  Destroy Local Scratch Vector:

  CALL VecDestroy(FieldsScratch, Error)

  IF(U % MPI_Rank == 0)THEN
    PRINT*
    PRINT*, "  INFO: Finished Evolution"
    PRINT*, "    Cycles Evolved = ", iCycle
    PRINT*, "    Time           = ", Time
  END IF

  !  Compare Evolved Solution with Initial Condition in the L1 Norm:

  CALL DMDAVecGetArrayF90( &
         F % FieldsDA, InitialCondition, FieldsOld, Error)
  CALL DMDAVecGetArrayF90( &
         F % FieldsDA, F % FieldsLocal,  FieldsNew, Error)

  MyL1Norm = Zero
  DO iX3 = iBX3, iEX3
    DO iX2 = iBX2, iEX2
      DO iX1 = iBX1, iEX1

        MyL1Norm &
          = MyL1Norm &
              + ABS(FieldsNew(D, iX1, iX2, iX3) - FieldsOld(D, iX1, iX2, iX3))

      END DO
    END DO
  END DO

  CALL MPI_REDUCE( &
         MyL1Norm, L1Norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
         U % Communicator, Error)

  IF(U % MPI_Rank == 0)THEN
    PRINT*
    PRINT*, "  INFO: L1 Norm = ", L1Norm
  END IF

  CALL DMDAVecRestoreArrayF90( &
         F % FieldsDA, InitialCondition, FieldsOld, Error)
  CALL DMDAVecRestoreArrayF90( &
         F % FieldsDA, F % FieldsGlobal, FieldsNew, Error)

  CALL VecDestroy(InitialCondition, Error)

  !  Release Access to Coordinates:

  CALL DMDAVecRestoreArrayF90( &
         M % CoordinateDA, M % Positions % InnerEdgeGlobal, InnerEdge, Error)
  CALL DMDAVecRestoreArrayF90( &
         M % CoordinateDA, M % Positions % CenterGlobal,    Center,    Error)
  CALL DMDAVecRestoreArrayF90( &
         M % CoordinateDA, M % Positions % OuterEdgeGlobal, OuterEdge, Error)

  END ASSOCIATE ! iBX1, iBX2, iBX3, etc.

  NULLIFY(F)
  NULLIFY(M)

  CALL DestroyUniverse(U)

END PROGRAM Advection3D
