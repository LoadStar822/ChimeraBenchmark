from .base import Registry
from ..tools.bracken import BrackenTool
from ..tools.chimera import ChimeraTool
from ..tools.ganon import GanonTool
from ..tools.kraken2 import Kraken2Tool
from ..tools.sylph import SylphTool
from ..tools.taxor import TaxorTool

TOOLS = Registry()
TOOLS.register("bracken", BrackenTool)
TOOLS.register("chimera", ChimeraTool)
TOOLS.register("ganon", GanonTool)
TOOLS.register("kraken2", Kraken2Tool)
TOOLS.register("sylph", SylphTool)
TOOLS.register("taxor", TaxorTool)
