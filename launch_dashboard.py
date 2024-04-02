from contextlib import redirect_stdout
from io import StringIO
import webview
from pymaldiproc.dashboard import app
from pymaldiproc import VERSION

# TODO: Figure out how to launch this script from a batch file.
if __name__ == '__main__':
    stream = StringIO()
    with redirect_stdout(stream):
        window = webview.create_window(f'pyMALDIproc Dashboard {VERSION}', app.server)
        webview.start()
