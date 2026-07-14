CHANGELOG
=========

This update is the June AROMA package with the lift locked, the validation
extended, and the notes moved into per-component PDFs. The controller and the
direct multi-step predictor approach are unchanged. The plant model has one
correction, listed below.

Since the June package
----------------------
- Lift locked at n_z = 47, state-only. The source supply command was removed
  from the dictionary and now enters only as an input, which addresses the
  midterm remark about inputs in the dictionary.
- Plant correction. shared/network/physics/components/pipe_rhs_bidir.m: the
  advection speed in the reverse-flow branch of the upwind scheme used the
  signed velocity where it needs the magnitude (-(v/dx) -> -(abs(v)/dx)). The
  heating-mode runs in this package never enter that branch, so no shipped
  result changes, but the branch is now correct if it is ever exercised.
- The iterated one-step baseline reads delivered heat off the lift the same way
  the proposed controller does, so the comparison uses a strong baseline; the
  two Koopman controllers then tie at the flow-limited ceiling.
- Documentation. Each component folder carries a *_validation.pdf: what it does,
  the main equations, how to reproduce it, every figure discussed, and the coded
  checks with the assert threshold next to the observed value. README.md is the
  entry point; these replace the earlier per-component notes.
- Added. The four-controller benchmark (visualizations/fair_comparison,
  see comparison_validation.pdf), an independent QP verification
  (tests/validate_qp.m), a source-pulse transport demo (demo_pulse.m), the
  identification-data adequacy checks, and the dictionary study behind the
  47-feature choice.
- The five validators were tightened and now run in one pass (run_all_tests.m,
  5/5 PASS). Result files and figures were regenerated under the above; a few
  comment and consistency corrections carry no change to any shipped result.
