import os
import tempfile

class TempSettings:
  def __init__(self):
    self.setup()

  def setup(self):
    repo = tempfile.mkdtemp()
    #test, if the cmr-settings-file is available or not:
    os.environ["CMR_SETTINGS_FILE"] = os.path.join(repo, "cmr-settings-file")
    f = open(os.environ["CMR_SETTINGS_FILE"], "w")
    f.write("[user]\nauto_add_to_db=False\nplugin_paths=\ndefault_privacy=640\nlocation=cmr-test\n")
    f.close()

  def __del__(self):
    os.system("rm -f %s"%os.environ["CMR_SETTINGS_FILE"] )
    os.system("rmdir  %s"%os.path.dirname(os.environ["CMR_SETTINGS_FILE"]))
