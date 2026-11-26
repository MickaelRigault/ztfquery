"""Offline regression tests for https://github.com/MickaelRigault/ztfquery/issues/33

IRSA passwords containing special characters (&, =, #, +, %, spaces, unicode)
must be percent-encoded when building the JOSSO login URL, otherwise the
login request is malformed and authentication fails.

All HTTP traffic is monkeypatched: these tests never touch the network.
"""
from urllib.parse import parse_qs, urlparse

import pytest

from ztfquery import io


class _FakeResponse:
    cookies = {}


class _FakeSession:
    """Just enough of a requests.Session for get_cookie()."""

    def __init__(self):
        self.cookies = {}
        self.urls = []

    def get(self, url, **kwargs):
        self.urls.append(url)
        return _FakeResponse()


@pytest.fixture
def login_urls(monkeypatch):
    """Capture every URL that get_cookie() would send to the IRSA login service."""
    calls = []

    def fake_get(url, *args, **kwargs):
        calls.append(url)
        return _FakeResponse()

    monkeypatch.setattr(io.requests, "get", fake_get)
    return calls


@pytest.mark.parametrize(
    "username, password",
    [
        ("plainuser", "plainpass"),
        ("u&xer=1", "pa&ss=word"),
        ("plainuser", "100%!"),
        ("plainuser", "with space"),
        ("plainuser", "p@ss#word+1"),
        ("plainuser", "a%26b"),  # literal percent must be double-encoded
        ("plainuser", "café ☕"),  # unicode
        ("plainuser", "&"),  # worst case: only separators
    ],
)
def test_credentials_survive_url_roundtrip(login_urls, username, password):
    session = _FakeSession()
    io.get_cookie(username=username, password=password, session=session)

    # get_cookie() hits the login URL once through the module-level requests.get
    # (and once through the session when one is provided).
    assert login_urls, "expected an HTTP call to the login url"
    url = login_urls[0]

    parsed = urlparse(url)
    query = parse_qs(parsed.query)

    # Exactly the three expected parameters: a raw '&' in the credentials must
    # not split the query string into extra parameters.
    assert set(query) == {"josso_cmd", "josso_username", "josso_password"}

    # Decoding the query must recover the exact original credentials.
    assert query["josso_cmd"] == ["login"]
    assert query["josso_username"] == [username]
    assert query["josso_password"] == [password]


def test_login_url_path_and_host_unchanged(login_urls):
    io.get_cookie(username="u", password="p", session=_FakeSession())
    parsed = urlparse(login_urls[0])
    assert parsed.scheme == "https"
    assert parsed.netloc == "irsa.ipac.caltech.edu"
    assert parsed.path == "/account/signon/login.do"
