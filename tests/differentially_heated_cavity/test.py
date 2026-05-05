#!/usr/bin/env python
def test_differentially_heated_cavity():
    import numpy as np
    from read_single_field_binary import read_single_field_binary

    scal, xp, yp, zp, xu, yv, zw = read_single_field_binary(
        "sca_001_fld_0000020_b_001.bin",
        1,
        np.array([1, 1, 1]),
    )
    assert np.all(np.isfinite(scal))
    assert np.mean(scal[0, :, :]) < np.mean(scal[-1, :, :])

    vex, xp, yp, zp, xu, yv, zw = read_single_field_binary(
        "vex_fld_0000020_b_001.bin",
        1,
        np.array([1, 1, 1]),
    )
    vez, xp, yp, zp, xu, yv, zw = read_single_field_binary(
        "vez_fld_0000020_b_001.bin",
        1,
        np.array([1, 1, 1]),
    )
    assert np.all(np.isfinite(vex))
    assert np.all(np.isfinite(vez))


if __name__ == "__main__":
    test_differentially_heated_cavity()
    print("Passed!")
