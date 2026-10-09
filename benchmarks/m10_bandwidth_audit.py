#!/usr/bin/env python3
"""M10.3A: audit structural bandwidth of the transient Jacobian."""

from dataclasses import dataclass


@dataclass
class Pattern:
    n: int
    entries: set[tuple[int, int]]

    @property
    def size(self):
        return 2 * self.n + 2

    @property
    def lower(self):
        return max(i - j for i, j in self.entries)

    @property
    def upper(self):
        return max(j - i for i, j in self.entries)


def jacobian_pattern(n: int) -> Pattern:
    """Reproduce the index assignments in ASSEMBLE_JACOBIAN.

    Indices are one-based, matching the Fortran implementation.
    """
    if n < 2:
        raise ValueError("n must be >= 2")

    entries = {
        (1, 1),
        (2, 1),
        (2, 2),
        (2, 3),
    }

    for i in range(n):
        im = 2 * i + 3
        il = 1 if i == 0 else 2 * i + 2
        ir = 2 * i + 4
        row = 2 * i + 3

        # Momentum equation
        entries.update({
            (row, il),
            (row, im),
            (row, ir),
        })

        # Continuity equation
        entries.add((row + 1, ir))
        entries.add((row + 1, im))

        if i < n - 1:
            entries.add((row + 1, im + 2))

    size = 2 * n + 2
    assert all(
        1 <= row <= size and 1 <= col <= size
        for row, col in entries
    )

    return Pattern(n=n, entries=entries)


def print_pattern(pattern: Pattern):
    size = pattern.size

    print(f"\nJacobian pattern: N={pattern.n}, unknowns={size}")
    print("    " + "".join(str(j // 10 % 10) for j in range(1, size + 1)))
    print("    " + "".join(str(j % 10) for j in range(1, size + 1)))

    for i in range(1, size + 1):
        line = "".join(
            "X" if (i, j) in pattern.entries else "."
            for j in range(1, size + 1)
        )
        print(f"{i:3d} {line}")


def main():
    for n in (2, 3, 4, 10, 20, 40, 80, 100):
        pattern = jacobian_pattern(n)

        print(
            f"N={n:3d} "
            f"unknowns={pattern.size:3d} "
            f"nonzeros={len(pattern.entries):4d} "
            f"lower={pattern.lower} "
            f"upper={pattern.upper}"
        )

        assert pattern.lower == 2
        assert pattern.upper == 1

    print_pattern(jacobian_pattern(4))
    print("\nPASS M10.3A: structural Jacobian bandwidth")


if __name__ == "__main__":
    main()
