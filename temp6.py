from functools import singledispatch

@singledispatch
def save(x) -> None:
    raise TypeError(f"Unsupported type: {type(x)!r}")

@save.register
def _(x: int) -> None:
    print(f"int {x}")

@save.register
def _(x: float) -> None:
    print(f"float {x}")

# Usage
save(5)     # prints: 5
save(3.14)  # prints: float 3.14