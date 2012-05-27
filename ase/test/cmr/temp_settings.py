import os
import tempfile

class TempSettings:
    """Creates settings file for cmr"""
    def __init__(self):
        self.created = False
        self.setup()

    def setup(self):
        repo = tempfile.mkdtemp()
        #test, if the cmr-settings-file is available or not:
        os.environ["CMR_SETTINGS_FILE"] = os.path.join(repo, "cmr-settings-file")
        f = open(os.environ["CMR_SETTINGS_FILE"], "w")
        f.write("[user]\nauto_add_to_db=False\nplugin_paths=\ndefault_privacy=640\nlocation=cmr-test\n")
        f.close()
        self.created = True

    def __del__(self):
        if not self.created:
            return
        sf = os.environ["CMR_SETTINGS_FILE"]
        dr = os.path.dirname(sf)
        if os.path.exists(sf):
            os.unlink(os.environ["CMR_SETTINGS_FILE"])
        if os.path.exists(dr):
            os.rmdir(dr)
