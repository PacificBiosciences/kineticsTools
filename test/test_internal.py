import os
import pytest
import kineticsTools.internal.basic as B

def test_basic():
    expected = 'anything'
    assert expected == B.getIpdModelFilename(expected, 'foo', [])

    with pytest.raises(Exception) as excinfo:
        B.getIpdModelFilename(None, 'unknown', [])
    assert 'Chemistry cannot be identified' in str(excinfo.value)

    with pytest.raises(Exception) as excinfo:
        B.getIpdModelFilename(None, 'foo', [])
    assert 'No kinetics model available for this chemistry' in str(excinfo.value)


def test_path(monkeypatch):
    def isfile(fn):
        if fn in ('path1/foo.h5', 'path2/foo.h5'):
            return True
        if fn == 'pathmissing/foo.h5':
            return False
        raise Exception('Called! {!r}'.format(fn))
    monkeypatch.setattr(os.path, 'isfile', isfile)

    chem = 'foo'

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['pathmissing', 'path1'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1', 'pathmissing'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1', 'path2'])


def test_path_with_prefixed_chem(monkeypatch):
    def isfile(fn):
        if fn in ('path1/Sfoo.h5', 'path2/Sfoo.h5'):
            return True
        if fn == 'pathmissing/Sfoo.h5':
            return False
        raise Exception('Called! {!r}'.format(fn))
    monkeypatch.setattr(os.path, 'isfile', isfile)

    chem = 'S/foo' # S/ prefix is weird for now.

    expected = 'path1/Sfoo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1'])

    expected = 'path1/Sfoo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['pathmissing', 'path1'])

    expected = 'path1/Sfoo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1', 'pathmissing'])

    expected = 'path1/Sfoo.h5'
    assert expected == B.getIpdModelFilename(None, chem, ['path1', 'path2'])


def test_getResourcePathSpec(monkeypatch):
    monkeypatch.setenv('SMRT_CHEMISTRY_BUNDLE_DIR', 'foo')
    assert 'foo/kineticsTools:bar' == B.getResourcePathSpec('bar')
