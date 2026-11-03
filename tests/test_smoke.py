import fork


def test_package_imports() -> None:
    assert fork.__name__ == "fork"
