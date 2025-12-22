from chimera_bench.registry.base import Registry


def test_registry_register_and_get():
    reg = Registry()
    reg.register("x", 1)
    assert reg.get("x") == 1


def test_registry_duplicate_raises():
    reg = Registry()
    reg.register("x", 1)
    try:
        reg.register("x", 2)
        assert False, "expected ValueError"
    except ValueError:
        pass
