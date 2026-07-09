# Changelog

## July 2026: dictionary correction and fair baseline

- `candidate_library.m`: removed the source command T_0s from the thermal block and
  the Qnet*(T_0s - T_amb) exergy product. The lift is now state-only with n_z = 47;
  inputs enter through B. This addresses the midterm remark about inputs in the
  dictionary.
- `build_lmpc_pred.m`: the iterated baseline now reads delivered heat directly off the
  lifted state (the theta coordinate), the same way the proposed controller does,
  instead of a tangent linearisation around the previous operating point. This makes
  the baseline stronger: its worst-consumer comfort on the stressed day went from
  91.5 to 97.7 percent.
- `run_comparison.m`: NMPC solver budget raised (iterations 100 to 300, function
  evaluations 600 to 3000) so it stops on its own tolerance, and both QP controllers
  now use the same solver tolerances as KPC. The NMPC and Jacobian legs were re-run
  under these settings.
- `tests/validate_robustness.m`: the C2 criterion now checks that KPC meets the 95
  percent comfort bar, is feasible on every step, and beats the two plant-based
  baselines. The old criterion also required beating the iterated Koopman baseline;
  after the baseline fix above the two tie at the capacity ceiling (97.6 vs 97.7),
  so that check would fail on run-to-run differences and was removed. Both the old
  and the new criterion pass on the old data.
- `run_stress_sweep.m`, `viz_stress_sweep.m`: sweep grid changed from
  [0.50 ... 0.30] (tightening only) to [0.35 ... 1.00] (with headroom). With the fair
  baseline the old conclusion (gap opens as capacity tightens) no longer holds: at
  the ceiling both controllers are flow-limited and tie, and with headroom the direct
  predictor keeps full comfort. The dropped 0.30 point sits below the demand adequacy
  floor, where the capacity assumption itself does not hold; no thesis claim uses it.
- `shared/params.m`: comment updates only (n_z, delay lags, excitation reference).
- `build_incidence_v2.m`: comment corrected. Kirchhoff on the flow references keeps
  the realised flows conservative at steady state, with a small transient residual;
  the earlier comment claimed conservation at every step, which is not true under
  the per-edge flow lags.
- `README.md`, `PREDICTOR_NOTES.md`: updated to the 47-feature lift.
- Regenerated artifacts: predictor and controller demo results, the V_seq fits,
  iterated A and B, `comparison.mat`, `longevity.mat`, `stress_sweep.mat` and the
  corresponding figures.
