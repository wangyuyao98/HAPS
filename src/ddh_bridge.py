"""Python bridge used by R(reticulate) for Dynamic-DeepHit-style prediction.

Default backend is "external", which expects a user-provided python module:
    module: args["external_module"] (default "ddh_backend")
    function: args["external_function"] (default "fit_predict")

That external function should return survival predictions for calibration/test
on a common time grid. This bridge standardizes outputs for R.

An optional backend "km_debug" is included for smoke tests only.
"""

from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np


def _to_array_1d(x: Any, name: str) -> np.ndarray:
    arr = np.asarray(x).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} is empty.")
    return arr


def _to_matrix(x: Any, name: str) -> np.ndarray:
    arr = np.asarray(x)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array; got shape {arr.shape}.")
    if arr.shape[0] == 0 or arr.shape[1] == 0:
        raise ValueError(f"{name} is empty with shape {arr.shape}.")
    return arr


def _n_rows(dat: Any, id_col: str) -> int:
    """Number of rows in an R/pandas/dict-like table."""
    return int(np.asarray(dat[id_col]).reshape(-1).size)


def _import_external_backend(module_name: str, module_path: Optional[str] = None):
    """Import external backend module.

    First tries normal import; if that fails and module_name has no dot,
    tries loading from explicit module_path, then common local paths.
    """
    try:
        return importlib.import_module(module_name)
    except Exception:
        candidates = []
        if module_path:
            candidates.append(Path(module_path))
        # fallbacks robust to reticulate where __file__ may be undefined
        candidates.append(Path.cwd() / "python" / f"{module_name}.py")
        candidates.append(Path.cwd() / f"{module_name}.py")
        if "__file__" in globals():
            candidates.append(Path(__file__).resolve().with_name(f"{module_name}.py"))

        for candidate in candidates:
            try:
                cand = candidate.resolve()
            except Exception:
                cand = candidate
            if not cand.exists():
                continue
            spec = importlib.util.spec_from_file_location(module_name, str(cand))
            if spec is None or spec.loader is None:
                continue
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod

        raise


def _km_curve_rel(durations: np.ndarray, events: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Simple KM step curve in relative time (>= 0)."""
    durations = np.asarray(durations, dtype=float)
    events = np.asarray(events, dtype=int)
    valid = np.isfinite(durations) & (durations >= 0)
    durations = durations[valid]
    events = events[valid]
    if durations.size == 0:
        return np.array([0.0]), np.array([1.0])

    # event times only
    t_unique = np.sort(np.unique(durations[events == 1]))
    if t_unique.size == 0:
        return np.array([0.0]), np.array([1.0])

    surv_vals = []
    s = 1.0
    for t in t_unique:
        n_risk = np.sum(durations >= t)
        d_t = np.sum((durations == t) & (events == 1))
        if n_risk > 0:
            s *= (1.0 - d_t / n_risk)
        surv_vals.append(s)
    return t_unique, np.asarray(surv_vals, dtype=float)


def fit_predict_ddh(
    dat_tr: Any,
    dat_cal: Any,
    dat_test: Any,
    dat_tr_tau: Any,
    dat_cal_tau: Any,
    dat_test_tau: Any,
    tau: float,
    start_col: str,
    stop_col: str,
    event_col: str,
    id_col: str,
    covars: Optional[Iterable[str]] = None,
    args: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    args = dict(args or {})
    backend = str(args.get("backend", "external"))
    covars = list(covars or [])

    if backend == "km_debug":
        # Debug-only backend to verify R/Python plumbing.
        t_rel = np.asarray(dat_tr_tau[stop_col], dtype=float) - float(tau)
        e = np.asarray(dat_tr_tau[event_col], dtype=int)
        time_rel, surv_base = _km_curve_rel(t_rel, e)
        surv_tr = np.tile(surv_base, (_n_rows(dat_tr_tau, id_col), 1))
        surv_cal = np.tile(surv_base, (_n_rows(dat_cal_tau, id_col), 1))
        surv_test = np.tile(surv_base, (_n_rows(dat_test_tau, id_col), 1))
        return {
            "pred_time": float(tau) + time_rel,
            "pred_surv_tr": surv_tr,
            "pred_surv_cal": surv_cal,
            "pred_surv_test": surv_test,
            "id_tr": np.asarray(dat_tr_tau[id_col]),
            "id_cal": np.asarray(dat_cal_tau[id_col]),
            "id_test": np.asarray(dat_test_tau[id_col]),
            "backend": backend,
        }

    if backend != "external":
        raise ValueError(
            f"Unknown ddh backend '{backend}'. Use 'external' or 'km_debug'."
        )

    module_name = str(args.get("external_module", "ddh_backend"))
    module_path = args.get("external_module_path", None)
    fn_name = str(args.get("external_function", "fit_predict"))

    try:
        mod = _import_external_backend(module_name, module_path=module_path)
    except Exception as exc:
        raise ImportError(
            "Could not import external Dynamic-DeepHit backend module "
            f"'{module_name}'. Please provide args$external_module and ensure "
            "it is installed in the python environment used by reticulate."
        ) from exc

    if not hasattr(mod, fn_name):
        raise AttributeError(
            f"External module '{module_name}' has no function '{fn_name}'."
        )
    fn = getattr(mod, fn_name)

    res = fn(
        dat_tr=dat_tr,
        dat_cal=dat_cal,
        dat_test=dat_test,
        dat_tr_tau=dat_tr_tau,
        dat_cal_tau=dat_cal_tau,
        dat_test_tau=dat_test_tau,
        tau=float(tau),
        start_col=start_col,
        stop_col=stop_col,
        event_col=event_col,
        id_col=id_col,
        covars=covars,
        args=args,
    )
    if not isinstance(res, dict):
        raise ValueError("External ddh backend must return a dict.")

    pred_time = res.get("pred_time", None)
    if pred_time is None:
        time_grid_rel = res.get("time_grid_rel", None)
        if time_grid_rel is None:
            raise ValueError(
                "External ddh output must include either 'pred_time' or 'time_grid_rel'."
            )
        pred_time = float(tau) + _to_array_1d(time_grid_rel, "time_grid_rel")
    else:
        pred_time = _to_array_1d(pred_time, "pred_time")

    surv_tr = res.get("pred_surv_tr", res.get("surv_tr", None))
    surv_cal = res.get("pred_surv_cal", res.get("surv_cal", None))
    surv_test = res.get("pred_surv_test", res.get("surv_test", None))
    if surv_tr is None or surv_cal is None or surv_test is None:
        raise ValueError(
            "External ddh output must include training/calibration/test survival matrices "
            "under keys pred_surv_tr/pred_surv_cal/pred_surv_test "
            "(or surv_tr/surv_cal/surv_test)."
        )

    surv_tr = _to_matrix(surv_tr, "pred_surv_tr")
    surv_cal = _to_matrix(surv_cal, "pred_surv_cal")
    surv_test = _to_matrix(surv_test, "pred_surv_test")
    if (
        surv_tr.shape[1] != pred_time.size
        or surv_cal.shape[1] != pred_time.size
        or surv_test.shape[1] != pred_time.size
    ):
        raise ValueError(
            "Column mismatch: survival matrices must align with pred_time length."
        )

    return {
        "pred_time": pred_time,
        "pred_surv_tr": surv_tr,
        "pred_surv_cal": surv_cal,
        "pred_surv_test": surv_test,
        "id_tr": res.get("id_tr", np.asarray(dat_tr_tau[id_col])),
        "id_cal": res.get("id_cal", np.asarray(dat_cal_tau[id_col])),
        "id_test": res.get("id_test", np.asarray(dat_test_tau[id_col])),
        "backend": backend,
    }
