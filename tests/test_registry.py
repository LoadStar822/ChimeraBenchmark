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

from chimera_bench.registry import TOOLS
from chimera_bench.tools.chimera import ChimeraTool
from chimera_bench.tools.ganon import GanonTool
from chimera_bench.tools.taxor import TaxorTool


def test_registry_has_chimera():
    tool = TOOLS.get("chimera")
    assert tool is ChimeraTool


def test_registry_has_ganon():
    tool = TOOLS.get("ganon")
    assert tool is GanonTool


def test_registry_has_taxor():
    tool = TOOLS.get("taxor")
    assert tool is TaxorTool
