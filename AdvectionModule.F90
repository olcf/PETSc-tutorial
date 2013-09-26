MODULE AdvectionModule

  USE RealKindModule, ONLY: &
    Half
  USE MeshModule, ONLY: &
    X1, X2, X3
  USE FieldsModule, ONLY: &
    D, V1, V2, V3

  IMPLICIT NONE
  PRIVATE

#include "finclude/petsc.h90"

  PUBLIC :: &
    AdvectExchangeCells3D, &
    AdvectNonExchangeCells3D

CONTAINS


  SUBROUTINE AdvectExchangeCells3D( &
               FieldsNew, FieldsOld, InnerEdge, OuterEdge, TimeStep, &
               iBX1, iBX2, iBX3, iEX1, iEX2, iEX3, SW)

    PetscReal, DIMENSION(:,:,:,:), POINTER :: &
      FieldsNew, &
      FieldsOld, &
      InnerEdge, &
      OuterEdge
    PetscReal, INTENT(in) :: &
      TimeStep
    PetscInt, INTENT(in) :: &
      iBX1, iBX2, iBX3, &
      iEX1, iEX2, iEX3, SW

    PetscInt :: &
      X1B, X1E, &
      X2B, X2E, &
      X3B, X3E

    !  1:
    X1B = iBX1
    X1E = iBX1 + SW - 1
    X2B = iBX2
    X2E = iEX2
    X3B = iBX3
    X3E = iEX3
    CALL Advect3D( &
           FieldsNew(D,                        &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           FieldsOld(D,                        &
             X1B - SW : X1E + SW,              &
             X2B - SW : X2E + SW,              &
             X3B - SW : X3E + SW),             &
           FieldsNew(V1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           OuterEdge(X1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E)  &
           - InnerEdge(X1,                     &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           TimeStep,                           &
           X1B = X1B, X1E = X1E,               &
           X2B = X2B, X2E = X2E,               &
           X3B = X3B, X3E = X3E, SW = SW)

    !  2:
    X1B = iEX1 - SW + 1
    X1E = iEX1
    X2B = iBX2
    X2E = iEX2
    X3B = iBX3
    X3E = iEX3
    CALL Advect3D( &
           FieldsNew(D,                        &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           FieldsOld(D,                        &
             X1B - SW : X1E + SW,              &
             X2B - SW : X2E + SW,              &
             X3B - SW : X3E + SW),             &
           FieldsNew(V1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           OuterEdge(X1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E)  &
           - InnerEdge(X1,                     &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           TimeStep,                           &
           X1B = X1B, X1E = X1E,               &
           X2B = X2B, X2E = X2E,               &
           X3B = X3B, X3E = X3E, SW = SW)

    !  3:
    X1B = iBX1 + SW
    X1E = iEX1 - SW
    X2B = iBX2
    X2E = iBX2 + SW - 1
    X3B = iBX3
    X3E = iEX3
    CALL Advect3D( &
           FieldsNew(D,                        &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           FieldsOld(D,                        &
             X1B - SW : X1E + SW,              &
             X2B - SW : X2E + SW,              &
             X3B - SW : X3E + SW),             &
           FieldsNew(V1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           OuterEdge(X1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E)  &
           - InnerEdge(X1,                     &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           TimeStep,                           &
           X1B = X1B, X1E = X1E,               &
           X2B = X2B, X2E = X2E,               &
           X3B = X3B, X3E = X3E, SW = SW)

    !  4:
    IF(iEX2 > iBX2)THEN
      X1B = iBX1 + SW
      X1E = iEX1 - SW
      X2B = iEX2 - SW + 1
      X2E = iEX2
      X3B = iBX3
      X3E = iEX3
      CALL Advect3D( &
             FieldsNew(D,                        &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             FieldsOld(D,                        &
               X1B - SW : X1E + SW,              &
               X2B - SW : X2E + SW,              &
               X3B - SW : X3E + SW),             &
             FieldsNew(V1,                       &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             OuterEdge(X1,                       &
               X1B : X1E, X2B : X2E, X3B : X3E)  &
             - InnerEdge(X1,                     &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             TimeStep,                           &
             X1B = X1B, X1E = X1E,               &
             X2B = X2B, X2E = X2E,               &
             X3B = X3B, X3E = X3E, SW = SW)
    END IF

    !  5:
    X1B = iBX1 + SW
    X1E = iEX1 - SW
    X2B = iBX2 + SW
    X2E = iEX2 - SW
    X3B = iBX3
    X3E = iBX3 + SW - 1
    CALL Advect3D( &
           FieldsNew(D,                        &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           FieldsOld(D,                        &
             X1B - SW : X1E + SW,              &
             X2B - SW : X2E + SW,              &
             X3B - SW : X3E + SW),             &
           FieldsNew(V1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           OuterEdge(X1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E)  &
           - InnerEdge(X1,                     &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           TimeStep,                           &
           X1B = X1B, X1E = X1E,               &
           X2B = X2B, X2E = X2E,               &
           X3B = X3B, X3E = X3E, SW = SW)

    !  6:
    IF(iEX3 > iBX3)THEN
      X1B = iBX1 + SW
      X1E = iEX1 - SW
      X2B = iBX2 + SW
      X2E = iEX2 - SW
      X3B = iEX3 - SW + 1
      X3E = iEX3
      CALL Advect3D( &
             FieldsNew(D,                        &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             FieldsOld(D,                        &
               X1B - SW : X1E + SW,              &
               X2B - SW : X2E + SW,              &
               X3B - SW : X3E + SW),             &
             FieldsNew(V1,                       &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             OuterEdge(X1,                       &
               X1B : X1E, X2B : X2E, X3B : X3E)  &
             - InnerEdge(X1,                     &
               X1B : X1E, X2B : X2E, X3B : X3E), &
             TimeStep,                           &
             X1B = X1B, X1E = X1E,               &
             X2B = X2B, X2E = X2E,               &
             X3B = X3B, X3E = X3E, SW = SW)
    END IF

  END SUBROUTINE AdvectExchangeCells3D


  SUBROUTINE AdvectNonExchangeCells3D( &
               FieldsNew, FieldsOld, InnerEdge, OuterEdge, TimeStep, &
               iBX1, iBX2, iBX3, iEX1, iEX2, iEX3, SW)

    PetscReal, DIMENSION(:,:,:,:), POINTER :: &
      FieldsNew, &
      FieldsOld, &
      InnerEdge, &
      OuterEdge
    PetscReal, INTENT(in) :: &
      TimeStep
    PetscInt, INTENT(in) :: &
      iBX1, iBX2, iBX3, &
      iEX1, iEX2, iEX3, SW

    PetscInt :: &
      X1B, X1E, &
      X2B, X2E, &
      X3B, X3E

    X1B = iBX1 + SW
    X1E = iEX1 - SW
    X2B = iBX2 + SW
    X2E = iEX2 - SW
    X3B = iBX3 + SW
    X3E = iEX3 - SW
    CALL Advect3D( &
           FieldsNew(D,                        &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           FieldsOld(D,                        &
             X1B - SW : X1E + SW,              &
             X2B - SW : X2E + SW,              &
             X3B - SW : X3E + SW),             &
           FieldsNew(V1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           OuterEdge(X1,                       &
             X1B : X1E, X2B : X2E, X3B : X3E)  &
           - InnerEdge(X1,                     &
             X1B : X1E, X2B : X2E, X3B : X3E), &
           TimeStep,                           &
           X1B = X1B, X1E = X1E,               &
           X2B = X2B, X2E = X2E,               &
           X3B = X3B, X3E = X3E, SW = SW)

  END SUBROUTINE AdvectNonExchangeCells3D


  SUBROUTINE Advect3D(D_NEW, D_OLD, V1, dX1, dt, &
               X1B, X1E, X2B, X2E, X3B, X3E, SW)

    PetscInt, INTENT(in) :: &
      X1B, X1E, X2B, X2E, X3B, X3E, SW
    PetscReal, DIMENSION(X1B : X1E, X2B : X2E, X3B : X3E), INTENT(inout) :: &
      D_NEW
    PetscReal, DIMENSION(X1B - SW : X1E + SW, &
                         X2B - SW : X2E + SW, &
                         X3B - SW : X3E + SW), INTENT(in) :: &
      D_OLD
    PetscReal, DIMENSION(X1B : X1E, X2B : X2E, X3B : X3E), INTENT(in) :: &
      V1, &
      dX1
    PetscReal, INTENT(in) :: &
      dt

    PetscReal :: &
      FluxOuterX1, &
      FluxInnerX1
    PetscInt :: &
      iX1, iX2, iX3

    DO iX3 = X3B, X3E
      DO iX2 = X2B, X2E
        DO iX1 = X1B, X1E

          FluxInnerX1 &
            = Half * ( ( V1(iX1, iX2, iX3) &
                         + ABS(V1(iX1, iX2, iX3)) ) * D_OLD(iX1-1, iX2, iX3) &
                       + ( V1(iX1, iX2, iX3) &
                           - ABS(V1(iX1, iX2, iX3)) ) * D_OLD(iX1, iX2, iX3) )

          FluxOuterX1 &
            = Half * ( ( V1(iX1, iX2, iX3) &
                         + ABS(V1(iX1, iX2, iX3)) ) * D_OLD(iX1, iX2, iX3) &
                       + ( V1(iX1, iX2, iX3) &
                           - ABS(V1(iX1, iX2, iX3)) ) * D_OLD(iX1+1, iX2, iX3) )

          D_NEW(iX1, iX2, iX3) &
            = D_OLD(iX1, iX2, iX3) &
                - dt * ( FluxOuterX1 - FluxInnerX1 ) / dX1(iX1, iX2, iX3)

        END DO
      END DO
    END DO

  END SUBROUTINE Advect3D


END MODULE AdvectionModule
