"""Dynamic-DeepHit-style backend for model.pred='ddh' using PyTorch.

This module implements a practical dynamic survival model that consumes the
full covariate trajectory up to prediction time tau and outputs subject-level
survival curves on a common time grid.

Interface contract:
    fit_predict(...) -> dict with keys
        - time_grid_rel OR pred_time
        - pred_surv_cal (n_cal x m)
        - pred_surv_test (n_test x m)
        - optional id_cal / id_test
"""

from __future__ import annotations

import random
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset
except Exception as exc:  # pragma: no cover - handled at runtime
    torch = None
    nn = None
    DataLoader = None
    TensorDataset = None
    _TORCH_IMPORT_ERROR = exc
else:
    _TORCH_IMPORT_ERROR = None


def _col(dat: Any, name: str) -> np.ndarray:
    """Read a column from dict-like or pandas-like objects as numpy array."""
    try:
        return np.asarray(dat[name])
    except Exception as exc:
        raise KeyError(f"Cannot read column '{name}' from input object.") from exc


def _scalar(x: Any) -> Any:
    if isinstance(x, np.generic):
        return x.item()
    return x


def _build_sequences(
    dat_full: Any,
    dat_tau: Any,
    tau: float,
    start_col: str,
    stop_col: str,
    id_col: str,
    covars: Sequence[str],
    include_time_features: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build padded trajectory tensors for ids in dat_tau from dat_full.

    Returns:
        X: (n_subject, max_len, n_feat)
        lengths: (n_subject,)
        ids_tau: (n_subject,)
    """
    ids_full = _col(dat_full, id_col)
    start_full = _col(dat_full, start_col).astype(float)
    stop_full = _col(dat_full, stop_col).astype(float)
    cov_full = {c: _col(dat_full, c).astype(float) for c in covars}

    # Stable ordering by id then start for trajectory extraction.
    ord_idx = np.lexsort((start_full, ids_full))
    ids_full = ids_full[ord_idx]
    start_full = start_full[ord_idx]
    stop_full = stop_full[ord_idx]
    for c in covars:
        cov_full[c] = cov_full[c][ord_idx]

    # Row-index map per subject in full long data.
    id_to_rows: Dict[Any, List[int]] = {}
    for i, sid in enumerate(ids_full):
        key = _scalar(sid)
        if key not in id_to_rows:
            id_to_rows[key] = []
        id_to_rows[key].append(i)

    ids_tau = _col(dat_tau, id_col)
    seq_list: List[np.ndarray] = []
    lengths: List[int] = []

    # Optional simple time encodings. Turn this off for a cleaner baseline-only
    # static comparison where future interval widths should not enter prediction.
    n_time_feat = 2 if include_time_features else 0
    n_feat = len(covars) + n_time_feat
    if n_feat == 0:
        n_feat = 1
    for sid in ids_tau:
        key = _scalar(sid)
        rows = id_to_rows.get(key, [])
        use_rows = [r for r in rows if start_full[r] <= float(tau) + 1e-12]

        if len(use_rows) == 0:
            # Defensive fallback: a single zero row.
            seq = np.zeros((1, n_feat), dtype=float)
            seq[0, -1] = 1e-6
        else:
            feat = np.zeros((len(use_rows), n_feat), dtype=float)
            for j, r in enumerate(use_rows):
                for cc, c in enumerate(covars):
                    feat[j, cc] = cov_full[c][r]
                if include_time_features:
                    feat[j, len(covars)] = start_full[r]
                    observed_stop = min(stop_full[r], float(tau))
                    feat[j, len(covars) + 1] = max(observed_stop - start_full[r], 1e-6)
                elif len(covars) == 0:
                    feat[j, 0] = 1.0
            seq = feat

        seq_list.append(seq)
        lengths.append(seq.shape[0])

    max_len = max(lengths) if len(lengths) > 0 else 1
    X = np.zeros((len(seq_list), max_len, n_feat), dtype=np.float32)
    for i, seq in enumerate(seq_list):
        X[i, : seq.shape[0], :] = seq.astype(np.float32)

    return X, np.asarray(lengths, dtype=np.int64), ids_tau


def _standardize_from_train(
    X_train: np.ndarray,
    len_train: np.ndarray,
    X_other: Sequence[Tuple[np.ndarray, np.ndarray]],
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Standardize sequence features with train-only moments."""
    n, t, d = X_train.shape
    mask_train = (np.arange(t)[None, :] < len_train[:, None]).reshape(n, t, 1)
    flat = X_train[mask_train.repeat(d, axis=2)].reshape(-1, d)
    mu = flat.mean(axis=0) if flat.size > 0 else np.zeros(d, dtype=np.float32)
    sd = flat.std(axis=0) if flat.size > 0 else np.ones(d, dtype=np.float32)
    sd = np.where(sd < 1e-6, 1.0, sd)

    def _apply(X: np.ndarray, L: np.ndarray) -> np.ndarray:
        Xs = (X - mu.reshape(1, 1, -1)) / sd.reshape(1, 1, -1)
        mask = np.arange(X.shape[1])[None, :] < L[:, None]
        Xs[~mask, :] = 0.0
        return Xs.astype(np.float32)

    train_std = _apply(X_train, len_train)
    others_std = [_apply(x, l) for x, l in X_other]
    return train_std, others_std


def _make_time_grid(train_dur: np.ndarray, args: Dict[str, Any]) -> np.ndarray:
    grid_rel = np.asarray(args.get("time_grid_rel", []), dtype=float).reshape(-1)
    grid_rel = grid_rel[np.isfinite(grid_rel) & (grid_rel > 0)]
    if grid_rel.size > 0:
        return np.unique(np.sort(grid_rel))

    n_bins = int(args.get("num_time_bins", 40))
    n_bins = max(5, n_bins)

    q = np.linspace(0, 1, n_bins + 1)[1:]
    grid = np.quantile(train_dur, q)
    grid = np.unique(np.asarray(grid, dtype=float))
    grid = grid[np.isfinite(grid) & (grid > 0)]
    if grid.size == 0:
        g = float(np.nanmax(train_dur)) if np.any(np.isfinite(train_dur)) else 1.0
        grid = np.array([max(g, 1e-3)], dtype=float)
    return grid


def _duration_to_bin(dur: np.ndarray, grid: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(grid, dur, side="left")
    idx = np.clip(idx, 0, len(grid) - 1)
    return idx.astype(np.int64)


class _DynamicDeepHitNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        num_bins: int,
        num_layers: int = 1,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        do_gru = dropout if num_layers > 1 else 0.0
        self.gru = nn.GRU(
            input_size=input_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True,
            dropout=do_gru,
        )
        self.head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, num_bins),
        )

    def forward(self, x: torch.Tensor, lengths: torch.Tensor) -> torch.Tensor:
        packed = nn.utils.rnn.pack_padded_sequence(
            x, lengths.cpu(), batch_first=True, enforce_sorted=False
        )
        _, h = self.gru(packed)
        z = h[-1]
        logits = self.head(z)
        hazards = torch.sigmoid(logits).clamp(1e-6, 1 - 1e-6)
        return hazards


def _disc_surv_nll(
    hazards: torch.Tensor, bins: torch.Tensor, events: torch.Tensor
) -> torch.Tensor:
    eps = 1e-8
    log_h = torch.log(hazards + eps)
    log_1mh = torch.log(1.0 - hazards + eps)
    csum = torch.cumsum(log_1mh, dim=1)
    i = torch.arange(hazards.shape[0], device=hazards.device)
    surv_before = torch.where(
        bins > 0,
        csum[i, bins - 1],
        torch.zeros_like(events, dtype=hazards.dtype),
    )
    loglik_event = surv_before + log_h[i, bins]
    loglik_cens = csum[i, bins]
    loglik = torch.where(events > 0.5, loglik_event, loglik_cens)
    return -loglik.mean()


def _predict_surv_matrix(
    model: _DynamicDeepHitNet,
    X: np.ndarray,
    lengths: np.ndarray,
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    model.eval()
    ds = TensorDataset(
        torch.from_numpy(X).float(),
        torch.from_numpy(lengths).long(),
    )
    dl = DataLoader(ds, batch_size=batch_size, shuffle=False)
    out = []
    with torch.no_grad():
        for xb, lb in dl:
            xb = xb.to(device)
            lb = lb.to(device)
            hz = model(xb, lb)
            sv = torch.cumprod(1.0 - hz, dim=1)
            out.append(sv.detach().cpu().numpy())
    surv = np.vstack(out) if len(out) > 0 else np.zeros((0, 0), dtype=float)
    surv = np.clip(np.minimum.accumulate(surv, axis=1), 0.0, 1.0)
    return surv


def _km_fallback(
    dat_tr_tau: Any,
    dat_cal_tau: Any,
    dat_test_tau: Any,
    tau: float,
    stop_col: str,
    event_col: str,
    id_col: str,
) -> Dict[str, Any]:
    """Fallback for degenerate folds (e.g., no events)."""
    t_rel = _col(dat_tr_tau, stop_col).astype(float) - float(tau)
    e = _col(dat_tr_tau, event_col).astype(int)
    evt = np.sort(np.unique(t_rel[e == 1]))
    if evt.size == 0:
        evt = np.array([max(np.nanmax(t_rel), 1e-3)], dtype=float)

    s = 1.0
    surv = []
    for t in evt:
        n_risk = np.sum(t_rel >= t)
        d_t = np.sum((t_rel == t) & (e == 1))
        if n_risk > 0:
            s *= (1.0 - d_t / n_risk)
        surv.append(s)
    surv = np.asarray(surv, dtype=float)

    n_cal = len(_col(dat_cal_tau, id_col))
    n_test = len(_col(dat_test_tau, id_col))
    n_tr = len(_col(dat_tr_tau, id_col))
    return {
        "time_grid_rel": evt,
        "pred_surv_tr": np.tile(surv, (n_tr, 1)),
        "pred_surv_cal": np.tile(surv, (n_cal, 1)),
        "pred_surv_test": np.tile(surv, (n_test, 1)),
        "id_tr": _col(dat_tr_tau, id_col),
        "id_cal": _col(dat_cal_tau, id_col),
        "id_test": _col(dat_test_tau, id_col),
    }


def fit_predict(
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
    """Train dynamic DeepHit-style model and predict survival curves."""
    if torch is None:
        raise ImportError(
            "PyTorch is required for ddh_backend.fit_predict but is not installed."
        ) from _TORCH_IMPORT_ERROR

    args = dict(args or {})
    covars = list(covars or [])
    seed = int(args.get("seed", 123))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    # Build trajectory tensors from full long data for risk-set ids at tau.
    include_time_features = bool(args.get("include_time_features", True))
    X_tr, L_tr, id_tr = _build_sequences(
        dat_tr, dat_tr_tau, tau, start_col, stop_col, id_col, covars,
        include_time_features=include_time_features,
    )
    X_ca, L_ca, id_ca = _build_sequences(
        dat_cal, dat_cal_tau, tau, start_col, stop_col, id_col, covars,
        include_time_features=include_time_features,
    )
    X_te, L_te, id_te = _build_sequences(
        dat_test, dat_test_tau, tau, start_col, stop_col, id_col, covars,
        include_time_features=include_time_features,
    )

    y_tr = _col(dat_tr_tau, stop_col).astype(float) - float(tau)
    e_tr = _col(dat_tr_tau, event_col).astype(int)
    y_ca = _col(dat_cal_tau, stop_col).astype(float) - float(tau)
    y_te = _col(dat_test_tau, stop_col).astype(float) - float(tau)
    y_tr = np.maximum(y_tr, 1e-8)
    y_ca = np.maximum(y_ca, 1e-8)
    y_te = np.maximum(y_te, 1e-8)

    # Fallback for degenerate folds.
    if len(y_tr) < 8 or np.sum(e_tr == 1) == 0:
        return _km_fallback(
            dat_tr_tau=dat_tr_tau,
            dat_cal_tau=dat_cal_tau,
            dat_test_tau=dat_test_tau,
            tau=tau,
            stop_col=stop_col,
            event_col=event_col,
            id_col=id_col,
        )

    grid_rel = _make_time_grid(y_tr, args)
    b_tr = _duration_to_bin(y_tr, grid_rel)

    X_tr, (X_ca, X_te) = _standardize_from_train(X_tr, L_tr, [(X_ca, L_ca), (X_te, L_te)])

    hidden_dim = int(args.get("hidden_size", 64))
    num_layers = int(args.get("num_layers", 1))
    dropout = float(args.get("dropout", 0.1))
    lr = float(args.get("lr", 1e-3))
    weight_decay = float(args.get("weight_decay", 1e-5))
    epochs = int(args.get("epochs", 80))
    batch_size = int(args.get("batch_size", 64))
    batch_size = max(4, batch_size)

    use_gpu = bool(args.get("use_gpu", False)) and torch.cuda.is_available()
    device = torch.device("cuda" if use_gpu else "cpu")

    model = _DynamicDeepHitNet(
        input_dim=X_tr.shape[2],
        hidden_dim=hidden_dim,
        num_bins=len(grid_rel),
        num_layers=num_layers,
        dropout=dropout,
    ).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    ds = TensorDataset(
        torch.from_numpy(X_tr).float(),
        torch.from_numpy(L_tr).long(),
        torch.from_numpy(b_tr).long(),
        torch.from_numpy(e_tr.astype(np.float32)).float(),
    )
    dl = DataLoader(ds, batch_size=batch_size, shuffle=True)

    model.train()
    for _ in range(epochs):
        for xb, lb, bb, eb in dl:
            xb = xb.to(device)
            lb = lb.to(device)
            bb = bb.to(device)
            eb = eb.to(device)
            hz = model(xb, lb)
            loss = _disc_surv_nll(hz, bb, eb)
            opt.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            opt.step()

    surv_tr = _predict_surv_matrix(model, X_tr, L_tr, batch_size=batch_size, device=device)
    surv_cal = _predict_surv_matrix(model, X_ca, L_ca, batch_size=batch_size, device=device)
    surv_test = _predict_surv_matrix(model, X_te, L_te, batch_size=batch_size, device=device)

    return {
        "time_grid_rel": grid_rel.astype(float),
        "pred_surv_tr": surv_tr.astype(float),
        "pred_surv_cal": surv_cal.astype(float),
        "pred_surv_test": surv_test.astype(float),
        "id_tr": id_tr,
        "id_cal": id_ca,
        "id_test": id_te,
    }
