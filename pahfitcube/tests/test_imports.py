"""These tests check if the right packages are installed by importing a
few modules.

"""


def test_run_pahfit_cube():
    from pahfitcube import run_pahfit_cube

    assert run_pahfit_cube is not None
