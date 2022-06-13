"""
WARNING: this version only works on *some* systems.
The version in get_notebook_info seems to always work.
"""

from jupyter_server import serverapp
from jupyter_server.utils import url_path_join
from pathlib import Path
import re
import requests
from IPython import get_ipython
from pathlib import Path



def get_notebook_path() -> str:
    """
    Return the name and path of the JupyterLab notebook invoking this function.
    Example:  "experiments/life_1D/diffusion/reach_equilibrium_1.ipynb"

    See https://stackoverflow.com/a/69096754/5478830 and
    https://github.com/jupyter/jupyter_client/pull/656

    :return:
    """
    kernelIdRegex = re.compile(r"(?<=kernel-)[\w\d\-]+(?=\.json)")
    
    kernelId = kernelIdRegex.search(get_ipython().config["IPKernelApp"]["connection_file"])[0]

    for jupServ in serverapp.list_running_servers():
        for session in requests.get(url_path_join(jupServ["url"], "api/sessions"), params={"token": jupServ["token"]}).json():
            if kernelId == session["kernel"]["id"]:
                return str(session["notebook"]['path'])
                # Path(jupServ["root_dir"]) / str(session["notebook"]['path'])
                # The above alt version would return something like "D:\Docs\[...]\experiments\life_1D\diffusion\reach_equilibrium_1.ipynb"



def get_notebook_basename() -> str:
    """
    Return the name - WITHOUT extension - of the JupyterLab notebook invoking this function.
    Example:  "reach_equilibrium_1"

    :return:
    """
    full_name = get_notebook_path()
    return Path(full_name).stem
    # Alternatively:  (os.path.splitext("/path/to/some/file.txt.zip.asc")[0])
