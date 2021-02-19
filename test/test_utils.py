import pytest
from atatutils.utils import str2select_cond


@pytest.mark.parametrize('query, want', [
    ('key1=value', [('key1', '=', 'value')]),
    ('key=true', [('key', '=', 1)]),
    ('key=false', [('key', '=', 0)]),
    ('key=false,temp=2', [('key', '=', 0), ('temp', '=', 2)]),
    ('key=false,temp=something', [('key', '=', 0), ('temp', '=', 'something')]),
    ('key>1', [('key', '>', 1)]),
    ('key<1', [('key', '<', 1)]),
    ('key<=1', [('key', '<=', 1)]),
    ('key>=1', [('key', '>=', 1)]),
])
def test_str2select_cond(query, want):
    got = str2select_cond(query)

    assert len(got) == len(want)
    print(got, want)
    for g, w in zip(got, want):
        assert all(x == y for x, y in zip(g, w))
