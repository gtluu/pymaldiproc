from contextlib import redirect_stdout
from io import StringIO
import webview
from pymaldiviz.dashboard import app
from pymaldiproc import VERSION


def main():
    stream = StringIO()
    with redirect_stdout(stream):
        webview.settings['ALLOW_DOWNLOADS'] = True
        window = webview.create_window(f'pyMALDIproc Dashboard {VERSION}', app.server)
        webview.start()


if __name__ == '__main__':
    main()
