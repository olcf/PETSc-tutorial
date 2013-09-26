MODULE UniverseModule

  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh
  USE FieldsModule, ONLY: &
    FieldsType, &
    CreateFields, &
    DestroyFields

  IMPLICIT NONE
  PRIVATE

#include "finclude/petsc.h90"

  TYPE, PUBLIC :: UniverseType
    CHARACTER(80) :: &
      Name
    PetscInt :: &
      Communicator, &
      MPI_Rank, &
      N_MPI_Ranks
    PetscLogDouble :: &
      WallTimeInitialize, &
      WallTimeFinalize
    TYPE(MeshType), POINTER :: &
      Mesh
    TYPE(FieldsType), POINTER :: &
      Fields
  END TYPE UniverseType

  PUBLIC :: &
    CreateUniverse, &
    DestroyUniverse

CONTAINS


  SUBROUTINE CreateUniverse(U, Name, &
               nZones, InnerBoundaries, OuterBoundaries, BoundaryConditions)

    TYPE(UniverseType), POINTER :: &
      U
    CHARACTER(len=*), INTENT(in) :: &
      Name
    PetscInt, DIMENSION(:), INTENT(in) :: &
      nZones
    PetscReal, DIMENSION(:), INTENT(in) :: &
      InnerBoundaries, &
      OuterBoundaries
    DMDABoundaryType, DIMENSION(:), INTENT(in) :: &
      BoundaryConditions

    PetscInt :: &
      Error

    CALL PetscInitialize(PETSC_NULL_CHARACTER, Error)

    ALLOCATE(U)
    U % Communicator = MPI_COMM_WORLD
    CALL MPI_Comm_Rank(U % Communicator, U % MPI_Rank,    Error)
    CALL MPI_Comm_Size(U % Communicator, U % N_MPI_Ranks, Error)

    U % Name = Name

    CALL PetscGetTime(U % WallTimeInitialize, Error)

    IF(U % MPI_RANK == 0)THEN
      PRINT*
      PRINT*, "INFO: Creating Universe ", TRIM(U % Name)
      PRINT*
      PRINT*, "  INFO: Number of MPI Ranks = ", U % N_MPI_Ranks
    END IF

    ALLOCATE(U % Mesh)
    CALL CreateMesh(U % Mesh, U % Communicator, U % MPI_Rank, &
           nZones, InnerBoundaries, OuterBoundaries, BoundaryConditions)

    ALLOCATE(U % Fields)
    CALL CreateFields(U % Fields, U % Mesh, U % Communicator, U % MPI_Rank)

  END SUBROUTINE CreateUniverse


  SUBROUTINE DestroyUniverse(U)

    TYPE(UniverseType), POINTER :: &
      U

    PetscInt :: &
      Error

    IF(U % MPI_Rank == 0)THEN
      PRINT*
      PRINT*, "INFO: Destroying Universe ", TRIM(U % Name)
    END IF

    CALL DestroyFields(U % Fields)
    DEALLOCATE(U % Fields)

    CALL DestroyMesh(U % Mesh)
    DEALLOCATE(U % Mesh)

    CALL PetscGetTime(U % WallTimeFinalize, Error)

    IF(U % MPI_Rank == 0)THEN
      PRINT*
      PRINT*, "  INFO: ", U % N_MPI_Ranks, " MPI Ranks Finished in ", &
        U % WallTimeFinalize - U % WallTimeInitialize, " Seconds"
      PRINT*
    END IF

    CALL PetscFinalize(Error)

    DEALLOCATE(U)

  END SUBROUTINE DestroyUniverse


END MODULE UniverseModule
