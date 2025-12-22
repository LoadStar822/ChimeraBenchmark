from .base import Registry
from ..tools.chimera import ChimeraTool

TOOLS = Registry()
TOOLS.register("chimera", ChimeraTool({"bin": "Chimera"}))
