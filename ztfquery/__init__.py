
__version__ ="0.3.2"

# - First time you call the script, it'll ask for your login info
from .tools import has_encryption, encrypt
if not has_encryption(): encrypt()
