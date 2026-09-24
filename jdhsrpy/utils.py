import re

__all__ = ["only_numbers"]

def only_numbers(string):
    """Remove all non-numeric characters from a string."""
    return re.sub(r"\D", "", string)