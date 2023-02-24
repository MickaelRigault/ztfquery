import os
from pathlib import Path
from cryptography.fernet import Fernet

ztfquery_encrypted_path = Path.cwd() / "tests" / ".ztfquery_crypt"
ztfquery_unencrypted_path = Path.home() / ".ztfquery"

key = os.getenv("ZTFQUERY_KEY")
fernet = Fernet(key)

with open(ztfquery_encrypted_path, "rb") as enc_file:
    encrypted = enc_file.read()

decrypted = fernet.decrypt(encrypted)

with open(ztfquery_unencrypted_path, "wb") as dec_file:
    dec_file.write(decrypted)
