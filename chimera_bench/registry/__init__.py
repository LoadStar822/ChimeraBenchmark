from .base import Registry
from ..tools.chimera import ChimeraTool
from ..tools.ganon import GanonTool

TOOLS = Registry()
TOOLS.register("chimera", ChimeraTool)
TOOLS.register("ganon", GanonTool)
