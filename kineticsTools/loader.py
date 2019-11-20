"""Basic functional stuff, I guess.
"""
import logging
import os

LOG = logging.getLogger(__name__)


def getIpdModelFilename(ipdModel, majorityChem, paramsPath):
    """
    ipdModel: str
    majorityChem: str
    """
    # In order of precedence they are:
    # 1. Explicit path passed to --ipdModel
    # 2. In-order through each directory listed in --paramsPath

    if ipdModel:
        LOG.info("Using passed-in kinetics model: {!r}".format(ipdModel))
        return ipdModel

    if majorityChem == 'unknown':
        msg = "Chemistry cannot be identified---cannot perform kinetic analysis"
        LOG.error(msg)
        raise Exception(msg)

    # Route any pre-Kiwi / pre-SSB Sequel chemistries to Seabiscuit training
    if majorityChem.startswith("S/P1") or majorityChem.startswith("S/P2"):
        majorityChem = "SP2-C2"
    # Route any post-SSB Sequel chemistries to Kiwi training (for now)
    elif majorityChem.startswith("S/"):
        majorityChem = "SP3-C3"

    # '/' is not a valid character in a file, unescaped--remove it
    majorityChem = majorityChem.replace("/", "")

    # go through each paramsPath in-order, checking if the model exists there or no
    for paramsPath in paramsPath:
        ipdModel = os.path.join(paramsPath, majorityChem + ".npz.gz")
        if os.path.isfile(ipdModel):
            LOG.info(
                "Using chemistry-matched kinetics model: {!r}".format(ipdModel))
            return ipdModel

    msg = "No kinetics model available for this chemistry ({!r}) on paramsPath {!r}".format(
        ipdModel, paramsPath)
    LOG.error(msg)
    raise Exception(msg)


def getResourcePathSpec(default_dir):
    """Create list of [${SMRT_CHEMISTRY_BUNDLE_DIR}/kineticsTools, {default_dir}].
    Return colon-separated string.
    """
    pths = []
    smrtChemBundlePath = os.environ.get("SMRT_CHEMISTRY_BUNDLE_DIR", None)
    if smrtChemBundlePath:
        LOG.info("found SMRT_CHEMISTRY_BUNDLE_DIR, prepending to default paramsPath")
        pths.append(os.path.join(smrtChemBundlePath, "kineticsTools"))
    pths.append(default_dir)
    return ':'.join(pths)
