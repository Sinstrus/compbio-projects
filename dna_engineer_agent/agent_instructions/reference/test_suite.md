# Test Suite

> Part of: [Agent Instructions](../README.md)

Comprehensive tests in `scripts/tools/tests/`:

1. **test_frame_offset.py** — BUG-001 (frame offset calculation)
2. **test_uniqueness_counting.py** — BUG-003 (double-strand counting)
3. **test_silent_classification.py** — Checkpoint 8 (silent mutations)
4. **test_parent_child_verification.py** — Checkpoint 10 (parent-child comparison)
5. **test_construct_verifier.py** — Full ConstructVerifier pipeline + KB integration
6. **test_design_reasoner.py** — DesignReasoner decision framework
7. **conftest.py** — Pytest fixtures and test data

**Run tests before major releases:**
```bash
python3 -m pytest scripts/tools/tests/ -v
```
