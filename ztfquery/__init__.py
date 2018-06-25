
__version__ ="0.5.1"

# - First time you call the script, it'll ask for your login info
from .io import has_encryption, encrypt
if not has_encryption(): encrypt()
