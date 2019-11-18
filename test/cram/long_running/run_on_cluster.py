#!/usr/bin/env python

import subprocess
import os.path

dir_name = os.path.dirname(__file__)
for org_code in ["Bsub", "Cagg", "Hpyl", "Mjan"]:
    script_name = "test_kineticsTools_%s.sh" % org_code
    with open(script_name, "w") as f:
        f.write("#!/bin/sh\n")
        f.write("module load cram\n")
        f.write("time cram %s/detect_and_identify_%s.t\n" %
                (dir_name, org_code))
    subprocess.call(["chmod", "755", "%s" % script_name])
    subprocess.call(["qsub", "-cwd", "-pe", "smp", "12", "%s" % script_name])
