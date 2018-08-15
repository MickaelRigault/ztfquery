
__version__ = "1.0.1"

# - First time you call the script, it'll ask for your login info
from .io import has_encryption, set_account
if not has_encryption():
    set_account("irsa")
