"""Basic functional stuff, I guess.
"""
import logging


def getIpdModelFilename(self):
    # In order of precedence they are:
    # 1. Explicit path passed to --ipdModel
    # 2. In-order through each directory listed in --paramsPath

    if self.args.ipdModel:
        logging.info("Using passed-in kinetics model: %s" % self.args.ipdModel)
        return self.args.ipdModel

    majorityChem = ReferenceUtils.loadAlignmentChemistry(self.alignments)
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
        logging.error("Chemistry cannot be identified---cannot perform kinetic analysis")
        sys.exit(1)

    # go through each paramsPath in-order, checking if the model exists there or no
    for paramsPath in self.args.paramsPath:
        ipdModel = os.path.join(paramsPath, majorityChem + ".h5")
        if os.path.isfile(ipdModel):
            logging.info("Using chemistry-matched kinetics model: {!r}".format(ipdModel))
            return ipdModel

    logging.error("Aborting, no kinetics model available for this chemistry: %s" % ipdModel)
    sys.exit(1)
