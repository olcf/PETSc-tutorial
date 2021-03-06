MODULE FieldsModule

  USE MeshModule, ONLY: &
    MeshType

  IMPLICIT NONE
  PRIVATE

#include "finclude/petsc.h90"

  PetscInt, PUBLIC, PARAMETER :: &
    D        = 0, &
    V1       = 1, &
    V2       = 2, &
    V3       = 3, &
    N_FIELDS = 4

  TYPE, PUBLIC :: FieldsType
    PetscInt :: &
      DegreesOfFreedom
    Vec :: &
      FieldsGlobal, &
      FieldsLocal
    DM :: &
      FieldsDA
  END TYPE FieldsType

  PUBLIC :: &
    CreateFields, &
    DestroyFields

CONTAINS


  SUBROUTINE CreateFields(F, M, Comm, Rank)

    TYPE(FieldsType), POINTER :: &
      F
    TYPE(MeshType), POINTER :: &
      M
    PetscInt, INTENT(in) :: &
      Comm, &
      Rank

    PetscErrorCode :: &
      Error

    F % DegreesOfFreedom = N_FIELDS

    IF(Rank == 0)THEN
      PRINT*
      PRINT*, "  INFO: Creating Fields"
      PRINT*, "    Degrees of Freedom = ", F % DegreesOfFreedom
    END IF

    SELECT CASE (M % nDimensions)
      CASE (1)
        CALL DMDACreate1D( &
               Comm,                 &
               M % BoundaryType(1),  &
               M % nZones(1),        &
               F % DegreesOfFreedom, &
               M % StencilWidth,     &
               PETSC_NULL_INTEGER,   &
               F % FieldsDA,         &
               Error)

      CASE (2)
        CALL DMDACreate2D( &
               Comm,                 &
               M % BoundaryType(1),  &
               M % BoundaryType(2),  &
               M % StencilType,      &
               M % nZones(1),        &
               M % nZones(2),        &
               PETSC_DECIDE,         &
               PETSC_DECIDE,         &
               F % DegreesOfFreedom, &
               M % StencilWidth,     &
               PETSC_NULL_INTEGER,   &
               PETSC_NULL_INTEGER,   &
               F % FieldsDA,         &
               Error)

      CASE (3)
        CALL DMDACreate3D( &
               Comm,                 &
               M % BoundaryType(1),  &
               M % BoundaryType(2),  &
               M % BoundaryType(3),  &
               M % StencilType,      &
               M % nZones(1),        &
               M % nZones(2),        &
               M % nZones(3),        &
               PETSC_DECIDE,         &
               PETSC_DECIDE,         &
               PETSC_DECIDE,         &
               F % DegreesOfFreedom, &
               M % StencilWidth,     &
               PETSC_NULL_INTEGER,   &
               PETSC_NULL_INTEGER,   &
               PETSC_NULL_INTEGER,   &
               F % FieldsDA,         &
               Error)

    END SELECT

    CALL DMCreateGlobalVector( &
           F % FieldsDA, F % FieldsGlobal, Error)

    CALL DMCreateLocalVector( &
           F % FieldsDA, F % FieldsLocal,  Error)

  END SUBROUTINE CreateFields


  SUBROUTINE DestroyFields(F)

    TYPE(FieldsType), POINTER :: &
      F

    PetscErrorCode :: &
      Error

    CALL VecDestroy(F % FieldsGlobal, Error)

    CALL VecDestroy(F % FieldsLocal,  Error)

    CALL DMDestroy(F % FieldsDA, Error)

  END SUBROUTINE DestroyFields


END MODULE FieldsModule
