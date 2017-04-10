"""Basic functional stuff, I guess.
"""
import logging
import os

LOG = logging.getLogger(__name__)


#majorityChem = ReferenceUtils.loadAlignmentChemistry(self.alignments)
def getIpdModelFilename(self, ipdModel, majorityChem):
    """
    ipdModel: str
    majorityChem: str
    """
    # In order of precedence they are:
    # 1. Explicit path passed to --ipdModel
    # 2. In-order through each directory listed in --paramsPath

    if ipdModel:
        logging.info("Using passed-in kinetics model: {!r}".format(ipdModel))
        return ipdModel

    # Temporary solution for Sequel chemistries: we do not
    # have trained kinetics models in hand yet for Sequel
    # chemistries.  However we have observed that the P5-C3
    # training seems to yield fairly good results on Sequel
    # chemistries to date.  So for the moment, we will use
    # that model for Sequel data.
    if majorityChem.startswith("S/"):
        logging.info("No trained model available yet for Sequel chemistries; modeling as P5-C3")
        majorityChem = "P5-C3"
    if majorityChem == 'unknown':
        msg = "Chemistry cannot be identified---cannot perform kinetic analysis"
        logging.error(msg)
        raise Exception(msg)

    # go through each paramsPath in-order, checking if the model exists there or no
    for paramsPath in self.args.paramsPath:
        ipdModel = os.path.join(paramsPath, majorityChem + ".h5")
        if os.path.isfile(ipdModel):
            logging.info("Using chemistry-matched kinetics model: {!r}".format(ipdModel))
            return ipdModel

    msg = "Aborting, no kinetics model available for this chemistry: {!r}".format(ipdModel)
    logging.error(msg)
    raise Exception(msg)
