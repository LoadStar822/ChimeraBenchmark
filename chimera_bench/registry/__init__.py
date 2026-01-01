from .base import Registry
from ..tools.chimera import ChimeraTool
from ..tools.ganon import GanonTool
from ..tools.taxor import TaxorTool

TOOLS = Registry()
TOOLS.register("chimera", ChimeraTool)
TOOLS.register("ganon", GanonTool)
TOOLS.register("taxor", TaxorTool)
