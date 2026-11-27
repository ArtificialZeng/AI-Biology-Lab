#!/usr/bin/env python3
"""Unit tests for Addendum B Stage 5A.3 CPU-vs-MPS lock policy."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

AUDIT_DIR = Path(__file__).resolve().parent
if str(AUDIT_DIR) not in sys.path:
    sys.path.insert(0, str(AUDIT_DIR))

from run_addb_gate5a3_mps import decide_addb_locked_device  # noqa: E402


class TestAddBGate5A3LockPolicy(unittest.TestCase):
    def test_lock_decision_thresholds(self) -> None:
        self.assertEqual(decide_addb_locked_device(1e-6, 5e-5, True, None), ("mps", True))
        self.assertEqual(decide_addb_locked_device(1e-5, 1e-4, True, None), ("mps", True))
        self.assertEqual(decide_addb_locked_device(1.01e-5, 5e-5, True, None), ("cpu", False))
        self.assertEqual(decide_addb_locked_device(1e-6, 1.01e-4, True, None), ("cpu", False))
        self.assertEqual(decide_addb_locked_device(0.0, 0.0, False, None), ("cpu", False))
        self.assertEqual(decide_addb_locked_device(0.0, 0.0, True, "err"), ("cpu", False))
        self.assertEqual(decide_addb_locked_device(float("nan"), 0.0, True, None), ("cpu", False))


if __name__ == "__main__":
    unittest.main()
