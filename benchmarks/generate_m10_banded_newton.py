#!/usr/bin/env python3
"""Generate the direct-banded Newton backend from the dense reference."""

from pathlib import Path

source = Path("src/fortran77/newton_solver.f")
target = Path("src/fortran77/newton_solver_banded.f")

text = source.read_text()


def replace_once(old, new):
    global text
    count = text.count(old)
    if count != 1:
        raise RuntimeError(
            f"Expected one occurrence, found {count}:\n{old}"
        )
    text = text.replace(old, new, 1)


replace_once(
    "SUBROUTINE NEWTON_SOLVE(",
    "SUBROUTINE NEWTON_SOLVE_BANDED("
)

replace_once(
    """      INTEGER NU,LDJ,I,J,LS,LINFO,ITER
      PARAMETER (LDJ=202)""",
    """      INTEGER NU,I,LS,LINFO,ITER
      INTEGER KL,KU,LDAB
      PARAMETER (KL=2,KU=1,LDAB=6)"""
)

replace_once(
    """      DOUBLE PRECISION R(202),RT(202),JAC(LDJ,202),AC(LDJ,202)
      DOUBLE PRECISION RHS(202),DU(202),UT(202),WORK(202)
      SAVE JAC,AC""",
    """      DOUBLE PRECISION R(202),RT(202),AB(LDAB,202)
      DOUBLE PRECISION RHS(202),DU(202),UT(202),WORK(202)"""
)

replace_once(
    """         CALL ASSEMBLE_JACOBIAN(N,WORK,DX,DT,THETA,D,A,T,Z,RS,
     &        MU,EPS,JAC,LDJ)
         DO 30 J=1,NU
            DO 20 I=1,NU
               AC(I,J) = JAC(I,J)
   20       CONTINUE
   30    CONTINUE
         DO 40 I=1,NU
            RHS(I) = -R(I)
   40    CONTINUE
         CALL SOLVE_DENSE(NU,AC,LDJ,RHS,DU,LINFO)""",
    """         CALL ASSEMBLE_JACOBIAN_BANDED(N,WORK,DX,DT,THETA,
     &        D,A,T,Z,RS,MU,EPS,AB,LDAB)
         DO 40 I=1,NU
            RHS(I) = -R(I)
   40    CONTINUE
         CALL SOLVE_BANDED(NU,AB,LDAB,KL,KU,RHS,DU,LINFO)"""
)

# Keep the generated source compatible with standard fixed-form Fortran.
long_lines = [
    (number, line)
    for number, line in enumerate(text.splitlines(), 1)
    if len(line) > 72
    and line
    and line[0] not in "Cc*! "
]

if long_lines:
    raise RuntimeError(f"Fixed-form lines exceed 72 columns: {long_lines}")

target.write_text(text)

print(f"Generated {target}")
print("Direct banded assembly: enabled")
print("Dense Jacobian allocation: removed")
print("Dense Jacobian copying: removed")
print("Dense linear solve: replaced")
