from certified_gaussian_theta import evaluate


def test_near_boundary_smoke():
    r = evaluate("1e-6", a="0", D=100, guard=20)
    assert r.selected_regime == "transformed"
    assert r.planned_transformed_depth <= 1
    assert float(r.certified_digits) > 100


def test_supported_rational_shift():
    r = evaluate("1", a="1/3", D=100, guard=20)
    assert r.R_coef == "0.0"
    assert float(r.certified_digits) > 100
