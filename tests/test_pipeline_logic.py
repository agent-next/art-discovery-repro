"""Unit tests for the pure functions in the census pipeline Python steps.

pipeline/census/*.py are scripts (stems start with digits), loaded here via
importlib. All tests are offline and deterministic (fixed seeds).
"""

import importlib.util
import itertools
import random
import sys
from collections import Counter
from pathlib import Path

import pytest

CENSUS = Path(__file__).resolve().parents[1] / "pipeline" / "census"


def _load(filename: str):
    """Load a pipeline/*.py script as a module."""
    path = CENSUS / filename
    spec = importlib.util.spec_from_file_location(f"pipeline_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


f02 = _load("02_filter.py")
f04 = _load("04_classify.py")
f05 = _load("05_sample_neighborhoods.py")
f06 = _load("06_score_partners.py")

_ids = itertools.count()


def _uid(prefix: str) -> str:
    return f"{prefix}-{next(_ids)}"


# ---------------------------------------------------------------------------
# 02_filter.py
# ---------------------------------------------------------------------------


def _hit(**kw) -> "f02.Hit":
    base = {
        "target": "t",
        "tlen": 300,
        "query": "pfam_retron",
        "qlen": 250,
        "bitscore": 100.0,
        "coverage": 1.0,
        "rt_class": "retron",
    }
    base.update(kw)
    return f02.Hit(**base)


def test_parse_hmmsearch_domtbl_best_domain_per_target(tmp_path: Path):
    domtbl = tmp_path / "hits.domtblout"
    domtbl.write_text(
        "# HMMER 3.4 domtblout\n"
        "#\n"
        # t1: two domain lines; parser keeps the higher bitscore (205.0)
        "t1 - 300 pfam_retron - 250 1e-60 210.0 0.1 1 1 1e-60 1e-60 205.0 0.1 "
        "1 250 5 250 1 250 0.99 -\n"
        "t1 - 300 pfam_retron - 250 1e-20 80.0 0.1 2 2 1e-20 1e-20 80.0 0.1 "
        "1 125 5 130 1 125 0.90 -\n"
        # t2: hmm 1..100 over model length 200 -> profile coverage 0.5;
        # pfam_dgr -> DGR class (coverage uses MODEL columns, not env coords)
        "t2 - 200 pfam_dgr - 200 1e-30 100.0 0.1 1 1 1e-30 1e-30 99.0 0.1 "
        "1 100 1 200 1 200 0.95 -\n"
        "malformed line\n"
        "\n"
    )
    hits = {h.target: h for h in f02.parse_hmmsearch_domtbl(domtbl)}
    assert set(hits) == {"t1", "t2"}
    assert hits["t1"].bitscore == 205.0  # best domain line kept
    assert hits["t1"].coverage == pytest.approx(1.0)
    assert hits["t1"].rt_class == "retron"
    assert hits["t2"].bitscore == 99.0
    assert hits["t2"].coverage == pytest.approx(0.5)
    assert hits["t2"].rt_class == "DGR"


def test_passes_filters_core_coverage_boundary():
    # paper: retain hits covering >=0.75 of the RT core profile
    assert f02.passes_filters(_hit(coverage=0.75)) is True
    assert f02.passes_filters(_hit(coverage=0.74)) is False


def test_passes_filters_weak_hit_rules():
    # paper Methods p.29: weak = bitscore <25 AND coverage <0.35 AND no YxDD,
    # jointly (grok review 2026-09-24: an OR reading discards real candidates)
    weak = dict(bitscore=24.9, coverage=0.30, has_yxdd=False)
    assert f02.is_weak_hit(_hit(**weak), f02.Thresholds()) is True
    # any single condition alone is NOT weak
    assert f02.is_weak_hit(_hit(bitscore=24.9, coverage=1.0, has_yxdd=True),
                           f02.Thresholds()) is False
    assert f02.is_weak_hit(_hit(bitscore=90.0, coverage=0.10, has_yxdd=False),
                           f02.Thresholds()) is False
    assert f02.is_weak_hit(_hit(bitscore=90.0, coverage=1.0, has_yxdd=False),
                           f02.Thresholds()) is False
    # boundary: bitscore exactly 25 is not < 25
    assert f02.is_weak_hit(_hit(bitscore=25.0, coverage=0.30, has_yxdd=False),
                           f02.Thresholds()) is False
    # and weak hits do not survive the filter set
    assert f02.passes_filters(_hit(**weak)) is False
    assert f02.passes_filters(_hit(bitscore=24.9, coverage=1.0, has_yxdd=True)) is True


def test_passes_filters_class_min_length():
    # paper: class-specific minimum length 225-346 aa (placeholder floor 225)
    assert f02.passes_filters(_hit(tlen=225)) is True
    assert f02.passes_filters(_hit(tlen=224)) is False
    # unknown class falls back to the default minimum
    assert f02.passes_filters(_hit(tlen=224, rt_class="mystery")) is False
    assert f02.passes_filters(_hit(tlen=225, rt_class="mystery")) is True


def test_classify_profile_keywords():
    assert f02.classify_profile("Pfam_retron_RT") == "retron"
    assert f02.classify_profile("group_II_intron_RT") == "group II intron"
    assert f02.classify_profile("CRISPR_cas_RT") == "CRISPR-assoc"
    assert f02.classify_profile("RVT_1") == "unknown"


# ---------------------------------------------------------------------------
# 04_classify.py
# ---------------------------------------------------------------------------


def _t1(profile: str, rt_class: str, bitscore: float) -> "f04.Tier1Hit":
    return f04.Tier1Hit(target="t", profile=profile, rt_class=rt_class,
                        bitscore=bitscore)


def _t2(bitscore: float = 150.0, coverage: float = 0.5,
        identity: float = 0.30) -> "f04.Tier2Hit":
    return f04.Tier2Hit(target="t", label="retron", bitscore=bitscore,
                        coverage=coverage, identity=identity)


def test_tier1_assign_clear_margin():
    hits = [_t1("pA", "retron", 150.0), _t1("pB", "DGR", 130.0)]
    assert f04.tier1_assign(hits) == "retron"


def test_tier1_margin_boundary():
    # paper: >=10-bit margin over the next-best profile of a different class
    hits_at = [_t1("pA", "retron", 150.0), _t1("pB", "DGR", 140.0)]
    assert f04.tier1_assign(hits_at) == "retron"  # margin exactly 10 -> assign
    hits_under = [_t1("pA", "retron", 150.0), _t1("pB", "DGR", 140.1)]
    assert f04.tier1_assign(hits_under) is None  # margin <10 -> defer


def test_tier1_floor_and_same_class_runner_up():
    # best profile below its class bitscore floor -> defer
    assert f04.tier1_assign([_t1("pA", "retron", 99.9)]) is None
    assert f04.tier1_assign([_t1("pA", "retron", 100.0)]) == "retron"
    # a runner-up of the SAME class does not count toward the margin
    hits = [_t1("pA", "retron", 150.0), _t1("pA2", "retron", 149.0),
            _t1("pB", "DGR", 100.0)]
    assert f04.tier1_assign(hits) == "retron"
    assert f04.tier1_assign([]) is None


def test_tier1_custom_floors():
    floors = {"retron": 200.0}
    assert f04.tier1_assign([_t1("pA", "retron", 150.0)], floors=floors) is None


def test_tier2_accept_boundaries():
    # paper: accept at bitscore >=150, coverage >=0.5, identity >=30%
    assert f04.tier2_accept(_t2()) is True
    assert f04.tier2_accept(_t2(bitscore=149.9)) is False
    assert f04.tier2_accept(_t2(coverage=0.49)) is False
    assert f04.tier2_accept(_t2(identity=0.29)) is False


def test_tier3_assign_nearest_leaves_agree():
    # paper: assign when the 5 nearest labeled leaves agree
    assert f04.tier3_assign(["retron"] * 5) == "retron"
    assert f04.tier3_assign(["retron"] * 4 + ["DGR"]) is None
    assert f04.tier3_assign(["retron"] * 4) is None  # fewer than 5 leaves


def test_classify_target_tier_precedence():
    assert f04.classify_target("t", [_t1("pA", "retron", 150.0)]) == (
        "tier1", "retron")
    assert f04.classify_target("t", [], _t2()) == ("tier2", "retron")
    assert f04.classify_target(
        "t", [], _t2(bitscore=10.0), ["UG"] * 5) == ("tier3", "UG")
    assert f04.classify_target("t", []) == ("unplaced", "unplaced")


def test_expected_class_counts_match_paper_census():
    assert len(f04.EXPECTED_CLASS_COUNTS) == 9
    assert sum(f04.EXPECTED_CLASS_COUNTS.values()) == 198_290


# ---------------------------------------------------------------------------
# 05_sample_neighborhoods.py
# ---------------------------------------------------------------------------


def _anchors(rt_class: str, n: int, clade: str = "") -> list:
    return [f05.Anchor(rt_id=_uid(rt_class), rt_class=rt_class, clade=clade)
            for _ in range(n)]


def test_sample_exhaustive_returns_all():
    members = _anchors("Abi", 5)
    assert f05.sample_exhaustive(members) == members


def test_sample_random_caps_at_quota():
    members = _anchors("DGR", 10)
    rng = random.Random(7)
    assert len(f05.sample_random(members, 3, rng)) == 3
    assert len(f05.sample_random(members, 100, rng)) == 10  # caps at pool


def test_plan_anchors_exhaustive_and_random_quotas():
    pools = {
        "Abi": f05.ClassPool(members=_anchors("Abi", 5)),
        "group II-like": f05.ClassPool(members=_anchors("group II-like", 4)),
        "DGR": f05.ClassPool(members=_anchors("DGR", 6)),
        "group II intron": f05.ClassPool(members=_anchors("group II intron", 600)),
    }
    sel = f05.plan_anchors(pools, random.Random(1))
    assert len(sel["Abi"]) == 5             # paper: exhaustive
    assert len(sel["group II-like"]) == 4   # paper: exhaustive
    assert len(sel["DGR"]) == 6             # quota 1,000 over a pool of 6
    assert len(sel["group II intron"]) == 500  # paper: random quota caps at 500


def test_plan_anchors_ug_proportional_sums_to_700():
    clades = {
        "u1": _anchors("UG", 500, clade="u1"),
        "u2": _anchors("UG", 300, clade="u2"),
        "u3": _anchors("UG", 200, clade="u3"),
    }
    pools = {"UG": f05.ClassPool(clades=clades)}
    sel = f05.plan_anchors(pools, random.Random(1))
    assert len(sel["UG"]) == 700  # paper: UG 700 proportional to clade size
    by_clade = Counter(a.clade for a in sel["UG"])
    assert by_clade == {"u1": 350, "u2": 210, "u3": 140}


def test_plan_anchors_novel_per_clade():
    clades = {"n1": _anchors("novel", 100, clade="n1"),
              "n2": _anchors("novel", 30, clade="n2")}
    pools = {"novel": f05.ClassPool(clades=clades)}
    sel = f05.plan_anchors(pools, random.Random(1))
    assert len(sel["novel"]) == 60 + 30  # paper: 60/clade, capped at clade size


def test_plan_anchors_unplaced_and_retron_control():
    clades = {f"L{i}": _anchors("unplaced", i + 1, clade=f"L{i}") for i in range(3)}
    unplaced = f05.ClassPool(clades=clades, singletons=_anchors("unplaced", 5))
    retron = f05.ClassPool(
        members=_anchors("retron", 10),
        controls=[f05.Anchor(rt_id="ctrl-1", rt_class="retron", is_control=True)],
    )
    sel = f05.plan_anchors({"unplaced": unplaced, "retron": retron},
                           random.Random(1))
    # one anchor per lineage (3 < 440) + all 5 singletons (< 200)
    assert len(sel["unplaced"]) == 3 + 5
    # paper: retron 1,000 random + 1 control
    assert len(sel["retron"]) == 10 + 1
    assert sel["retron"][-1].is_control


def test_plan_anchors_deterministic_same_seed():
    pools = {
        "DGR": f05.ClassPool(members=_anchors("DGR", 100)),
        "group II intron": f05.ClassPool(members=_anchors("group II intron", 600)),
        "UG": f05.ClassPool(clades={"u1": _anchors("UG", 500, clade="u1"),
                                    "u2": _anchors("UG", 500, clade="u2")}),
    }
    first = f05.plan_anchors(pools, random.Random(123))
    second = f05.plan_anchors(pools, random.Random(123))

    def ids(sel):
        return {k: [a.rt_id for a in v] for k, v in sel.items()}

    assert ids(first) == ids(second)


def test_locus_window_clamps_and_extends():
    # paper: neighborhood = locus +/- 10 kb flanks
    loc = f05.Locus(source="ENA", contig="c1", start=5_000, end=6_000)
    assert f05.locus_window(loc) == (0, 16_000)  # start clamps at 0
    loc2 = f05.Locus(source="NCBI", contig="c2", start=20_000, end=21_000)
    assert f05.locus_window(loc2) == (10_000, 31_000)
    assert f05.locus_window(loc2, flank=500) == (19_500, 21_500)


# ---------------------------------------------------------------------------
# 06_score_partners.py
# ---------------------------------------------------------------------------


def _gene(**kw) -> "f06.NeighborGene":
    base = {
        "gene_id": _uid("g"),
        "locus_id": "l1",
        "family": "fam1",
        "rt_clade": "cladeA",
        "rt_class": "retron",
        "cluster_id_90": "c1",
        "biosample": "b1",
        "rt_adjacent": True,
        "same_strand_as_rt": True,
        "distance_to_adjacent_bp": 0,
    }
    base.update(kw)
    return f06.NeighborGene(**base)


def _background(n: int) -> list:
    """Genes that never qualify for RT proximity."""
    return [_gene(rt_adjacent=False, same_strand_as_rt=False,
                  distance_to_adjacent_bp=None) for _ in range(n)]


def test_assign_family_prefers_pfam_over_cluster():
    # paper: family assignment by best Pfam-A match, else MMseqs2 cluster
    assert f06.assign_family("PfamA", "clu9") == "PfamA"
    assert f06.assign_family(None, "clu9") == "cluster:clu9"
    assert f06.assign_family("", "clu9") == "cluster:clu9"


def test_conserved_in_clade_boundaries():
    # paper filter 1: >=10% of the clade's loci OR >=3 loci
    assert f06.conserved_in_clade(3, 100) is True    # >=3 loci rule
    assert f06.conserved_in_clade(2, 100) is False   # 2 loci and 2% < 10%
    assert f06.conserved_in_clade(2, 20) is True     # exactly 10%
    assert f06.conserved_in_clade(1, 10) is True     # exactly 10%
    assert f06.conserved_in_clade(1, 11) is False    # just under 10%
    assert f06.conserved_in_clade(3, 0) is True      # zero-total: loci rule
    assert f06.conserved_in_clade(2, 0) is False


def test_passes_conservation_any_one_clade_suffices():
    occs = [
        _gene(locus_id="a1", rt_clade="big"),
        _gene(locus_id="b1", rt_clade="small"),
        _gene(locus_id="b2", rt_clade="small"),
    ]
    totals = {"big": 1000, "small": 20}  # 2/20 = 10% only in the small clade
    assert f06.passes_conservation(occs, totals) is True
    # distinct-locus counting: 3 occurrences on 2 loci of a large clade fail
    occs2 = [_gene(locus_id="l1"), _gene(locus_id="l1"), _gene(locus_id="l2")]
    assert f06.passes_conservation(occs2, {"cladeA": 100}) is False
    # ...but pass at >=3 distinct loci
    occs3 = [_gene(locus_id=f"l{i}") for i in range(3)]
    assert f06.passes_conservation(occs3, {"cladeA": 100}) is True


def test_passes_independence_each_boundary():
    # paper filter 2: >=3 clusters (90% id), >=3 biosamples, >=2 RT classes
    base = [
        _gene(cluster_id_90="c1", biosample="b1", rt_class="retron"),
        _gene(cluster_id_90="c2", biosample="b2", rt_class="DGR"),
        _gene(cluster_id_90="c3", biosample="b3", rt_class="retron"),
    ]
    assert f06.passes_independence(base) is True
    two_clusters = [base[0], base[1],
                    _gene(cluster_id_90="c1", biosample="b3")]
    assert f06.passes_independence(two_clusters) is False   # 2 < 3 clusters
    two_bios = [base[0], base[1],
                _gene(cluster_id_90="c3", biosample="b1")]
    assert f06.passes_independence(two_bios) is False       # 2 < 3 biosamples
    one_class = [_gene(cluster_id_90="c1", biosample="b1", rt_class="retron"),
                 _gene(cluster_id_90="c2", biosample="b2", rt_class="retron"),
                 _gene(cluster_id_90="c3", biosample="b3", rt_class="retron")]
    assert f06.passes_independence(one_class) is False      # 1 < 2 RT classes


def test_proximity_qualifies_rules():
    # paper filter 3 geometry: RT-adjacent AND (same strand OR <=100 bp away)
    assert f06.proximity_qualifies(_gene()) is True
    assert f06.proximity_qualifies(_gene(rt_adjacent=False)) is False
    assert f06.proximity_qualifies(
        _gene(same_strand_as_rt=False, distance_to_adjacent_bp=100)) is True
    assert f06.proximity_qualifies(
        _gene(same_strand_as_rt=False, distance_to_adjacent_bp=101)) is False
    assert f06.proximity_qualifies(
        _gene(same_strand_as_rt=False, distance_to_adjacent_bp=None)) is False


def test_permutation_pvalue_bounds_and_enrichment():
    occurrences = [_gene() for _ in range(5)]  # all qualify for proximity
    p = f06.permutation_pvalue(occurrences, _background(50),
                               n_permutations=200, rng=random.Random(1))
    assert 0.0 <= p <= 1.0
    # strongly enriched vs an all-non-adjacent background -> small p
    assert p <= f06.PERMUTATION_P_MAX


def test_permutation_pvalue_unenriched_is_one():
    occurrences = _background(4)  # observed qualifying count = 0
    p = f06.permutation_pvalue(occurrences, _background(50),
                               n_permutations=100, rng=random.Random(1))
    assert p == 1.0  # every permutation ties or exceeds


def test_promotes_null_draws_from_same_loci():
    # Paper p.30 (verbatim): "genes of the family had to lie nearer to the RT
    # than randomly drawn genes OF THE SAME LOCI (P <= 0.05)" — the null pool
    # is the family's own loci, not the global neighborhood pool. When the
    # same-loci genes are equally RT-proximal, the family must NOT be
    # promoted, even though a global pool (mostly distant, non-qualifying
    # loci) would give a small p. FAILING-FIRST against the global-pool code.
    occurrences = [
        _gene(locus_id="l1", cluster_id_90="c1", biosample="b1",
              rt_class="retron"),
        _gene(locus_id="l2", cluster_id_90="c2", biosample="b2",
              rt_class="DGR"),
        _gene(locus_id="l3", cluster_id_90="c3", biosample="b3",
              rt_class="retron"),
    ]
    same_loci_qualifying = [
        _gene(locus_id=f"l{i}", rt_adjacent=True, same_strand_as_rt=True,
              distance_to_adjacent_bp=0)
        for i in (1, 2, 3) for _ in range(10)
    ]
    other_loci_nonqualifying = [
        _gene(locus_id=f"far{i}", rt_adjacent=False, same_strand_as_rt=False,
              distance_to_adjacent_bp=None)
        for i in range(500)
    ]
    background = same_loci_qualifying + other_loci_nonqualifying
    promoted, p = f06.promotes(occurrences, {"cladeA": 100}, background,
                               n_permutations=200, rng=random.Random(1))
    assert p == 1.0, p  # same-loci draws qualify as often as the family
    assert promoted is False


def test_promotes_all_three_filters():
    occurrences = [
        _gene(locus_id="l1", cluster_id_90="c1", biosample="b1",
              rt_class="retron"),
        _gene(locus_id="l2", cluster_id_90="c2", biosample="b2",
              rt_class="DGR"),
        _gene(locus_id="l3", cluster_id_90="c3", biosample="b3",
              rt_class="retron"),
    ]
    promoted, p = f06.promotes(occurrences, {"cladeA": 100}, _background(50),
                               n_permutations=200, rng=random.Random(1))
    assert promoted is True
    assert p <= f06.PERMUTATION_P_MAX


def test_promotes_fails_early_returns_p1():
    # fails filter 1 (single locus in a large clade) -> short-circuit p=1.0
    promoted, p = f06.promotes([_gene(locus_id="l1")], {"cladeA": 1000}, [])
    assert promoted is False and p == 1.0
    # conserved but not independent (one cluster) -> still (False, 1.0)
    occs = [_gene(locus_id=f"l{i}", cluster_id_90="c1") for i in range(3)]
    promoted, p = f06.promotes(occs, {"cladeA": 100}, [])
    assert promoted is False and p == 1.0


def test_novel_eighth_class_vs_scattered_unplaced():
    # Paper p.29 (verbatim): "A further 25,737 clusters formed a loosely
    # organized group on the tree apart from any labeled member and were
    # treated as an eighth class of 'novel' RTs. The 137,385 scattered
    # clusters that could not be placed formed the ninth class of 'unplaced'
    # RTs." classify_target never returned 'novel' (raw-fidelity audit).
    # INTERPRETED criterion: a majority of unlabeled neighbors in the k
    # nearest leaves = grouped apart from labeled members -> novel; labeled
    # but disagreeing neighbors -> scattered -> unplaced.
    t1, t2 = [], None  # no tier-1/2 signal
    # neighborhood dominated by other unplaced clusters -> novel
    out = f04.classify_target("x", t1, t2, ["unplaced"] * 3 + ["retron", "retron"])
    assert out == ("tier3", "novel"), out
    # labeled neighbors that disagree -> scattered unplaced
    out = f04.classify_target("y", t1, t2, ["retron", "DGR", "retron", "UG", "Abi"])
    assert out == ("unplaced", "unplaced"), out
    # five agreeing labeled neighbors still place (unchanged behavior)
    out = f04.classify_target("z", t1, t2, ["retron"] * 5)
    assert out == ("tier3", "retron"), out
    # no neighbor data at all stays unplaced (no evidence of a group)
    out = f04.classify_target("w", t1, t2, None)
    assert out == ("unplaced", "unplaced"), out
