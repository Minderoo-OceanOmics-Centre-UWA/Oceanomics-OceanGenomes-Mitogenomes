#!/usr/bin/env python3
"""Per-row savepoints for the LCA upload scripts.

A plain try/except around cursor.execute() looks like it isolates a bad row, but
it does not: PostgreSQL aborts the whole transaction on any error, so every
later statement fails with "current transaction is aborted" and the eventual
commit is downgraded to a rollback. The scripts still printed a per-row tally
and exited 0, so a single unrepresentable value silently cost a sample its
entire upload -- which is how 87 assemblies across batch-12 .. batch-20 ended up
with no LCA rows in the database while their tasks reported success.

Wrapping each row in a savepoint makes the per-row handling real: a failure
rolls back that row alone and the rest of the batch still commits.
"""

from contextlib import contextmanager

_SAVEPOINT = "row_guard"


@contextmanager
def row_savepoint(cursor):
    """Run one row's statement so that its failure costs only that row.

    Re-raises, so the caller keeps its own except block: this only guarantees
    the transaction is still usable by the time that block runs.
    """
    cursor.execute(f"SAVEPOINT {_SAVEPOINT}")
    try:
        yield
    except Exception:
        cursor.execute(f"ROLLBACK TO SAVEPOINT {_SAVEPOINT}")
        raise
    else:
        cursor.execute(f"RELEASE SAVEPOINT {_SAVEPOINT}")
