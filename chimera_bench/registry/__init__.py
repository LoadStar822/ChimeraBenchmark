from .base import Registry
from ..tools.chimera import ChimeraTool
from ..tools.ganon import GanonTool
from ..tools.sylph import SylphTool
from ..tools.taxor import TaxorTool

TOOLS = Registry()
TOOLS.register("chimera", ChimeraTool)
TOOLS.register("ganon", GanonTool)
TOOLS.register("sylph", SylphTool)
TOOLS.register("taxor", TaxorTool)
