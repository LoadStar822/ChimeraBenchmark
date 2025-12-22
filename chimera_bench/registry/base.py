class Registry:
    def __init__(self) -> None:
        self._items = {}

    def register(self, name: str, item) -> None:
        if name in self._items:
            raise ValueError(f"Duplicate registration: {name}")
        self._items[name] = item

    def get(self, name: str):
        return self._items[name]
