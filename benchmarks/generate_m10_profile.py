#!/usr/bin/env python3
"""Generate instrumented copies of M9E solver and M10.1 benchmark driver.

Production sources are read-only inputs. Generated files go under build/.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'src/fortran77/newton_solver.f'
DRIVER = ROOT / 'benchmarks/m10_cpu_baseline.f90'
OUT = ROOT / 'build'


def instrument_calls(source, name, slot, expected):
    # Match fixed-form CALL statements and their continuation lines.
    pattern = re.compile(
        r'(?m)^([ ]{6,}CALL ' + re.escape(name) +
        r'\([^\n]*(?:\n[ ]{5}&[^\n]*)*)'
    )
    def replace(match):
        return (f'      CALL PROFILE_START({slot})\n' + match.group(1)
                + f'\n      CALL PROFILE_STOP({slot})')
    modified, count = pattern.subn(replace, source)
    if count != expected:
        raise RuntimeError(f'{name}: expected {expected} calls, found {count}')
    return modified


def main():
    source = SOURCE.read_text()
    if source.count('      IMPLICIT NONE') != 1:
        raise RuntimeError('Unexpected NEWTON_SOLVE declaration')
    source = source.replace('      IMPLICIT NONE',
        '      USE M10_PROFILE, ONLY: PROFILE_START, PROFILE_STOP\n'
        '      IMPLICIT NONE', 1)

    for name, slot, expected in (
        ('ASSEMBLE_RESIDUAL', 1, 2),
        ('ASSEMBLE_JACOBIAN', 2, 1),
        ('SOLVE_DENSE', 4, 1),
    ):
        source = instrument_calls(source, name, slot, expected)

    copy_loop = ("         DO 30 J=1,NU\n"
                 "            DO 20 I=1,NU\n"
                 "               AC(I,J) = JAC(I,J)\n"
                 "   20       CONTINUE\n"
                 "   30    CONTINUE")
    if source.count(copy_loop) != 1:
        raise RuntimeError('Cannot locate Jacobian copy loop')
    source = source.replace(copy_loop,
        '         CALL PROFILE_START(3)\n' + copy_loop +
        '\n         CALL PROFILE_STOP(3)', 1)

    driver = DRIVER.read_text()
    if driver.count('  implicit none') != 1:
        raise RuntimeError('Unexpected M10.1 driver declaration')
    driver = driver.replace('  implicit none',
        '  use m10_profile, only: profile_report\n  implicit none', 1)
    marker = '  elapsed = end_cpu-start_cpu'
    if driver.count(marker) != 1:
        raise RuntimeError('Cannot locate benchmark elapsed time')
    driver = driver.replace(marker, marker + '\n  call profile_report(elapsed)', 1)

    OUT.mkdir(exist_ok=True)
    (OUT / 'm10_newton_profile.f').write_text(source)
    (OUT / 'm10_cpu_profile.f90').write_text(driver)
    print('Generated M10.2 instrumented sources in build/')


if __name__ == '__main__':
    main()
