from .base import Registry
from ..tools.bracken import BrackenTool
from ..tools.chimera import ChimeraTool
from ..tools.centrifuger import CentrifugerTool
from ..tools.ganon import GanonTool
from ..tools.kraken2 import Kraken2Tool
from ..tools.metabuli import MetabuliTool
from ..tools.metaphlan import MetaPhlAnTool
from ..tools.sylph import SylphTool

TOOLS = Registry()
TOOLS.register("bracken", BrackenTool)
TOOLS.register("centrifuger", CentrifugerTool)
TOOLS.register("chimera", ChimeraTool)
TOOLS.register("ganon", GanonTool)
TOOLS.register("kraken2", Kraken2Tool)
TOOLS.register("metabuli", MetabuliTool)
TOOLS.register("metaphlan", MetaPhlAnTool)
TOOLS.register("sylph", SylphTool)
