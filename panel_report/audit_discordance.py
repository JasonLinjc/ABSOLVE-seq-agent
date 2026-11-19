"""Reconcile discordant manuscript calls and recover identifiable binomial counts."""
import argparse
import csv
import gzip
import subprocess
from tempfile import TemporaryDirectory
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from openpyxl import load_workbook
from scipy.optimize import least_squares
from scipy.special import expit

UPSTREAM_COMMIT = "271355dfe40cf46b952e82b88dd4521f23f9feb4"
FOLDERS = {"Unfiltered": "raw", "Dedup": "raw_dedup", "Dedud": "dedud_filtered",
           "Dedud+Dedup": "dedud_filtered_dedup"}
DONORS = ["D3734", "D3742", "D3908"]
ARMS = ["NoEP", "HiFi", "WT"]
TAGS = [d + "-" + ("WTCas9" if a == "WT" else a) for d in DONORS for a in ARMS]
X = np.array([[1, a == "HiFi", a == "WT", d == "D3742", d == "D3908"]
              for d in DONORS for a in ARMS], dtype=float)
STRATEGIES = ["Unfiltered", "Dedup", "Dedud+Dedup", "Dedud"]
TABLES = dict(zip(STRATEGIES, [3, 4, 5, 2]))


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, rows):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def bh(p):
    p = np.asarray(p)
    idx = np.argsort(p)
    q = np.minimum.accumulate((p[idx] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.minimum(1, q)[np.argsort(idx)]


def recover(f):
    # Invert the nine link-scale prediction variances for the nine cell weights.
    p = np.array([f["predicted-" + k] for k in TAGS])
    var = np.array([f["stdError-" + k] ** 2 for k in TAGS])

    def covariance(t):
        w = np.exp(t)
        return X @ np.linalg.solve(X.T @ (w[:, None] * X), X.T)

    def residual(t):
        return np.log(covariance(t).diagonal() / var)

    def jac(t):
        h = covariance(t)
        return -(h ** 2) * np.exp(t)[None, :] / h.diagonal()[:, None]

    opt = least_squares(residual, np.log(0.5 / var), jac=jac,
                        ftol=1e-13, xtol=1e-13, gtol=1e-13)
    assert opt.success and np.max(np.abs(opt.fun)) < 1e-9 and np.linalg.matrix_rank(opt.jac) == 9
    w = np.exp(opt.x)
    raw_n = w / (p * (1 - p))
    n = np.rint(raw_n).astype(int)
    assert max(abs(n - raw_n)) < 0.01 and min(n) > 0
    # Adjusted-score identities identify arm and donor event margins, not cell events.
    adjusted = (n * p - w * covariance(opt.x).diagonal() * (0.5 - p)).reshape(3, 3)
    event_arm = np.rint(adjusted.sum(axis=0)).astype(int)
    event_donor = np.rint(adjusted.sum(axis=1)).astype(int)
    err = max(abs(adjusted.sum(axis=0) - event_arm).max(),
              abs(adjusted.sum(axis=1) - event_donor).max())
    assert err < 0.01 and sum(event_arm) == sum(event_donor)
    stats = np.r_[event_arm.sum(), event_arm[1:], event_donor[1:]]

    # Forward fit from integer sufficient statistics, starting away from the solution.
    def score(beta):
        mu = expit(X @ beta)
        weights = n * mu * (1 - mu)
        h = weights * np.diag(X @ np.linalg.solve(X.T @ (weights[:, None] * X), X.T))
        return stats - X.T @ (n * mu - h * (0.5 - mu))

    beta0 = np.array([np.log((stats[0] + 0.5) / (n.sum() - stats[0] + 0.5)), 0, 0, 0, 0])
    fit = least_squares(score, beta0, ftol=1e-13, xtol=1e-13, gtol=1e-13)
    assert fit.success and max(abs(fit.fun)) < 1e-6
    labels = ["(Intercept)", "treatHiFi", "treatWTCas9", "donorD3742", "donorD3908"]
    beta_error = max(abs(fit.x - [f["estim-" + l] for l in labels]))
    assert beta_error < 1e-4
    return n.reshape(3, 3), event_arm, event_donor, beta_error, err


def grouped(rows, field):
    return np.array([[sum(int(r[field]) for r in rows if r["arm"] == a and r["donor"] == d)
                      for a in ARMS] for d in DONORS])


def public_counts(upstream, target, folder):
    target_name = target.replace("_", "-")
    paths = sorted((upstream / "test/absolveseq_edits" / folder).glob(target_name + "_*.gz"))
    rows = []
    for path in paths:
        _, donor, label, replicate, _ = path.name.split("_", 4)
        assert label in ("NoEP", "HiFiCas9")
        edited = total = 0
        with gzip.open(path, "rt") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                count = int(row["#Reads"])
                assert row["Unedited"] in ("True", "False")
                total += count
                if row["Unedited"] == "False":
                    edited += count
        rows.append(dict(donor=donor, replicate=replicate, arm=label.replace("Cas9", ""),
                         indel=edited, denominator=total))
    return rows


def audit(run, workbooks, upstream, output):
    # The quick start uses a detached checkout. This also works in the analysis image without Git.
    head = upstream / ".git/HEAD"
    revision = head.read_text().strip() if head.is_file() else ""
    if not revision or revision.startswith("ref:"):
        revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if revision != UPSTREAM_COMMIT:
        raise ValueError(f"Expected upstream commit {UPSTREAM_COMMIT}, found {revision}")
    rows = read(run / "panel_comparison.tsv")
    counts = read(run / "library_counts.tsv")
    fields = defaultdict(dict)
    for r in read(run / "model_field_comparison.tsv"):
        fields[r["target"], r["strategy"]][r["field"]] = float(r["agent"])
    hashes = json.loads((run / "provenance.json").read_text())["submitted_workbooks"]
    paper = {}
    for s, num in TABLES.items():
        path = workbooks / f"SupplementaryTable{num}.xlsx"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == hashes[path.name]
        wb = load_workbook(path, read_only=True, data_only=True)
        assert len(wb.worksheets) == 1
        vals = iter(wb.active.values)
        next(vals)
        header = next(vals)
        records = [dict(zip(header, v)) for v in vals if v and v[0]]
        assert len(records) == 79
        paper.update({(v["barcode"], s): v for v in records})
        wb.close()

    discordant = []
    for s in STRATEGIES:
        for editor in ("WT", "HiFi"):
            group = [r for r in rows if r["strategy"] == s and r["editor"] == editor]
            assert len(group) == len({r["target"] for r in group}) == 79
            arm = "WTCas9" if editor == "WT" else editor
            p_pub = np.array([paper[r["target"], s]["asymPval-treat" + arm] for r in group])
            p_agent = np.array([float(r["p_one_sided"]) for r in group])
            qa, qm = bh(p_agent), bh(p_pub)
            for i, r in enumerate(group):
                f = paper[r["target"], s]
                dm = 100 * (f["predicted-" + arm] - f["predicted-NoEP"])
                da = float(r["delta_pp"])
                np.testing.assert_allclose([dm, qm[i], qa[i]], [float(r[k]) for k in
                    ("manuscript_delta_pp", "manuscript_q_bh", "q_bh")], rtol=1e-10, atol=1e-11)
                yes_a, yes_m = da > 0.1 and qa[i] < 0.05, dm > 0.1 and qm[i] < 0.05
                assert (yes_a == yes_m) == (r["hit_agrees"] == "TRUE")
                if yes_a == yes_m:
                    continue
                swap_a, swap_m = p_agent.copy(), p_pub.copy()
                swap_a[i], swap_m[i] = p_pub[i], p_agent[i]
                delta_flip, q_flip = (dm > .1) != (da > .1), (qm[i] < .05) != (qa[i] < .05)
                discordant.append(dict(target=r["target"], strategy=s, editor=editor, table=TABLES[s],
                    hit_in="agent" if yes_a else "manuscript", gate=("delta+q" if q_flip else "delta") if delta_flip else "q",
                    agent_delta_pp=da, manuscript_delta_pp=dm, agent_q=qa[i], manuscript_q=qm[i],
                    agent_p=p_agent[i], manuscript_p=p_pub[i],
                    q_agent_p_other_78_manuscript=bh(swap_m)[i], q_manuscript_p_other_78_agent=bh(swap_a)[i],
                    agent_control_pct=r["control_rate_pct"], manuscript_control_pct=100*f["predicted-NoEP"],
                    agent_edited_pct=r["edited_rate_pct"], manuscript_edited_pct=100*f["predicted-"+arm]))
    discordant.sort(key=lambda r: (0 if r["editor"] == "WT" and r["strategy"] in STRATEGIES[:2] else 1,
                                   STRATEGIES.index(r["strategy"]), r["editor"], r["target"]))
    if not discordant:
        print("All 632 calls agree.")
        return
    cells, libraries, qc = [], [], []
    for r in discordant:
        t, s, ed = r["target"], r["strategy"], r["editor"]
        local = [v for v in counts if v["target"] == t and v["strategy"] == s]
        archived = public_counts(upstream, t, FOLDERS[s])
        key = lambda row: (row["donor"], row["replicate"], row["arm"])
        assert len(local) == len(set(map(key, local))) == 26
        archived_by_key = {key(v): v for v in archived}
        assert len(archived) == len(archived_by_key) == 17
        assert archived_by_key.keys() == {key(v) for v in local if v["arm"] != "WT"}
        values = {k: grouped(local, k) for k in ("indel", "denominator")}
        if s in STRATEGIES[:2]:
            assert all(v[k] == v["pre_step6_"+k] for v in local for k in values)
        an, ay, ad, ae, ar = recover(fields[t,s])
        np.testing.assert_array_equal(an, values["denominator"])
        np.testing.assert_array_equal(ay, values["indel"].sum(axis=0))
        np.testing.assert_array_equal(ad, values["indel"].sum(axis=1))
        mn, my, md, me, mr = recover(paper[t,s])
        pub_n, pub_y = (grouped(archived, k) for k in ("denominator","indel"))
        public_matches = (np.array_equal(pub_n[:, :2], mn[:, :2])
                          and np.array_equal(pub_y.sum(0)[:2], my[:2]))
        before_match = ""
        if s == "Dedud+Dedup":
            raw = public_counts(upstream, t, "dedud_dedup")
            assert len(raw) == 17
            before_match = bool(np.array_equal(grouped(raw, "denominator")[:, :2], mn[:, :2]) and
                                np.array_equal(grouped(raw, "indel").sum(0)[:2], my[:2]))
        # Only attribute cell event counts when the external control/HiFi table agrees with margins.
        cell_y = pub_y.copy()
        cell_y[:,2] = md - pub_y[:,:2].sum(1)
        if public_matches:
            assert np.all((0 <= cell_y) & (cell_y <= mn)) and np.array_equal(cell_y.sum(0), my)
        for j,a in enumerate(ARMS):
            if a in ("NoEP",ed):
                tag = "control" if a == "NoEP" else "edited"
                r.update({f"agent_{tag}_indel":ay[j], f"agent_{tag}_total":an[:,j].sum(),
                          f"manuscript_inferred_{tag}_indel":my[j],
                          f"manuscript_inferred_{tag}_total":mn[:,j].sum()})
            for i,d in enumerate(DONORS):
                cells.append(dict(target=t,strategy=s,editor=ed,donor=d,arm=a,
                    agent_indel=values["indel"][i,j],agent_total=an[i,j],
                    manuscript_inferred_total=mn[i,j], public_indel=pub_y[i,j] if j<2 else "",
                    public_total=pub_n[i,j] if j<2 else "",
                    manuscript_compatible_indel=cell_y[i,j] if public_matches else "",
                    event_source=("public allele table" if j<2 else "donor event margin minus public NoEP and HiFi") if public_matches else "unassigned"))
        r["public_control_hifi_match"] = public_matches
        qc.append(dict(target=t,strategy=s,agent_forward_beta_max_error=ae,manuscript_forward_beta_max_error=me,
                       agent_event_integer_max_error=ar,manuscript_event_integer_max_error=mr,
                       public_model_margins_match=public_matches,public_before_step6_margins_match=before_match))
        for v in local:
            p = archived_by_key.get(key(v), {})
            libraries.append({**v,"editor":ed,"public_indel":p.get("indel",""),"public_denominator":p.get("denominator","")})
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temp:
        stage = Path(temp) / "audit"
        stage.mkdir()
        for name, records in (("discordant_calls", discordant), ("donor_counts", cells),
                              ("library_counts_check", libraries), ("reconstruction_qc", qc)):
            write(stage / (name + ".tsv"), records)
        stage.rename(output)
    for group, count in sorted(Counter((r["strategy"], r["editor"], r["gate"]) for r in discordant).items()):
        print(*group, count, sep="\t")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path, help="Completed 26-library manuscript run directory")
    parser.add_argument("--workbooks", required=True, type=Path, help="Directory containing submitted Tables 2–5")
    parser.add_argument("--upstream", required=True, type=Path, help="Pinned upstream ABSOLVE-seq checkout")
    parser.add_argument("--output", required=True, type=Path, help="New directory for the audit tables")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must be a new directory")
    audit(args.run, args.workbooks, args.upstream, args.output)


if __name__ == "__main__":
    main()
