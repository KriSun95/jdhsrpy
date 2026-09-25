import re

__all__ = ["only_numbers"]

def only_numbers(string):
    """Remove all non-numeric characters from a string.
    
    For example, \"1a2H 3;\" becomes \"123\". 
    """
    return re.sub(r"\D", "", string)