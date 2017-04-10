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

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, 'foo', ['path1'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, 'foo', ['pathmissing', 'path1'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, 'foo', ['path1', 'pathmissing'])

    expected = 'path1/foo.h5'
    assert expected == B.getIpdModelFilename(None, 'foo', ['path1', 'path2'])

def test_getResourcePathSpec(monkeypatch):
    monkeypatch.setenv('SMRT_CHEMISTRY_BUNDLE_DIR', 'foo')
    assert 'foo/kineticsTools:bar' == B.getResourcePathSpec('bar')
