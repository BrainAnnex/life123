"""
When run in a Jupyterlab 3 notebook, it returns the notebook filename or the full path to the notebook.

Adapted from:  https://github.com/msm1089/ipynbname/blob/master/ipynbname/__init__.py

Added comments, and simplified to be more friendly to jupyterlab 3 - a simplification that
in some environments very dramatically speeds up execution: 
see https://github.com/msm1089/ipynbname/issues/2
and https://github.com/msm1089/ipynbname/pull/11

It may work with other versions of jupyterlab, or with the classic jupyter notebooks -
but not tested with them.
"""

import json
import urllib.error
import urllib.request
from pathlib import Path, PurePath
from typing import Generator, Tuple, Union

import ipykernel
from jupyter_core.paths import jupyter_runtime_dir
from traitlets.config import MultipleInstanceError


FILE_ERROR = "Can't identify the notebook {}."
CONN_ERROR = "Unable to access server;\n" \
           + "ipynbname requires either no security or token based security."



def _list_maybe_running_servers(pattern: str) -> Generator[dict, None, None]:
    """
    Iterate over the server info files of running notebook servers.

    :param pattern: Pattern to locate the names of interest.
                    On jupyterlab 3 on Windows, it is "jpserver-*.json";
                    on jupyterlab 3 on Linux, as well as on jupyterlab 2 and jupyter notebook, it appears to be "nbserver-*.json"
                    There seems to be inconsistency in the naming convention
    :return:
    """

    runtime_dir_str = jupyter_runtime_dir()
    # EXAMPLE on Windows: 'C:\\Users\\valerie\\AppData\\Roaming\\jupyter\\runtime'
    # EXAMPLE on Linux:   '/home/jovyan/.local/share/jupyter/runtime'
    
    runtime_dir = Path(runtime_dir_str)
    # EXAMPLE on Windows: WindowsPath('C:/Users/valerie/AppData/Roaming/jupyter/runtime')
    # EXAMPLE on Linux:   PosixPath('/home/jovyan/.local/share/jupyter/runtime')

    if runtime_dir.is_dir():
        # glob finds all the pathnames matching a specified pattern, and returns a generator object
        for file_name in runtime_dir.glob(pattern):
            # EXAMPLE of file_name on Windows:
            # 'C:\Users\valerie\AppData\Roaming\jupyter\runtime\jpserver-10368.json'
            # EXAMPLE on Linux:
            # '/home/jovyan/.local/share/jupyter/runtime/nbserver-25.json'
            yield json.loads(file_name.read_bytes())


            
def _get_kernel_id() -> str:
    """ Returns the kernel ID of the ipykernel.
    """
    connection_file_with_path = ipykernel.get_connection_file()
    # EXAMPLE on Windows:
    # 'C:\\Users\\valerie\\AppData\\Roaming\\jupyter\\runtime\\kernel-06c8322c-ef53-4d3e-a673-556a80bc0cc5.json'
    
    connection_file = Path(connection_file_with_path).stem
    kernel_id = connection_file.split('-', 1)[1]
    return kernel_id



def _get_sessions(srv):
    """ Given a server, returns sessions, or HTTPError if access is denied.
        NOTE: Works only when either there is no security or there is token
        based security. An HTTPError is raised if unable to connect to a 
        server.
    """
    try:
        qry_str = ""
        token = srv['token']
        if token:
            qry_str = f"?token={token}"
        url = f"{srv['url']}api/sessions{qry_str}"
        with urllib.request.urlopen(url) as req:
            return json.load(req)
    except Exception:
        raise urllib.error.HTTPError(CONN_ERROR)


        
def _find_nb_path(pattern: str) -> Union[Tuple[dict, PurePath], Tuple[None, None]]:
    """
    If successful in locating the notebook, return a pair consisting of:
    1) a dictionary of session data,
    2) a PurePath object for the located notebook.  EXAMPLE: PureWindowsPath('experiments/life_1D/diffusion/diffusion_1.ipynb')
    
    Otherwise, return (None, None)

    :param pattern: Pattern to locate the names of interest (for more info, see _list_maybe_running_servers)
    """
    try:
        kernel_id = _get_kernel_id()   # EXAMPLE:  'f191d796-178a-4a17-acf1-4c3933cd30fc'
    except (MultipleInstanceError, RuntimeError):
        return None, None    # Could not determine the kernel ID
    
    for srv in _list_maybe_running_servers(pattern):
        """EXAMPLE of srv: 
            {'base_url': '/', 
            'hostname': 'localhost', 
            'password': False, 
            'pid': 10368, 
            'port': 8888, 
            'root_dir': 'D:\\Docs\\- MY CODE\\BioSimulations\\life123-Win7', 
            'secure': False, 
            'sock': '', 
            'token': 'f0c1aa8a76721d96d27036f1f788d27c34eb9a5f661231e5', 
            'url': 'http://localhost:8888/', 
            'version': '1.17.0'}
        """
        try:
            sessions = _get_sessions(srv)
            for sess in sessions:
                """ EXAMPLE of sess:
                    {'id': '4c470050-021a-4000-b552-797e409e63a1', 
                    'path': 'experiments/life_1D/diffusion/diffusion_1.ipynb', 
                    'name': 'diffusion_1.ipynb', 
                    'type': 'notebook', 
                    'kernel': {'id': '4c8f0d11-2f99-46d4-8c89-e7198e10ac9f', 'name': 'python3', 'last_activity': '2022-06-13T03:31:45.436400Z', 'execution_state': 'idle', 'connections': 0}, 
                    'notebook': {'path': 'experiments/life_1D/diffusion/diffusion_1.ipynb', 'name': 'diffusion_1.ipynb'}}
                """
                if sess['kernel']['id'] == kernel_id:
                    return srv, PurePath(sess['notebook']['path'])
        except Exception:
            pass  # There may be stale entries in the runtime directory
        
    return None, None



def name() -> str:
    """
    Returns the short name of the notebook w/o the .ipynb extension,
    or raises a FileNotFoundError exception if it cannot be determined.
    """
    # First, look for names used on jupyterlab 3 on Windows (presumably, this is the new standard for names in jupyterlab)
    _, path = _find_nb_path("jpserver-*.json")
    if path:
        return path.stem

    # If nothing found thus far, look for names used on jupyterlab 3 on Linux, as well as on jupyterlab 2 and jupyter notebook
    _, path = _find_nb_path("nbserver-*.json")
    if path:
        return path.stem

    raise FileNotFoundError(FILE_ERROR.format('name'))


    
def path() -> Path:
    """ Returns the absolute path of the notebook,
        or raises a FileNotFoundError exception if it cannot be determined.
    """
    srv, path = _find_nb_path()
    if srv and path:
        root_dir = Path(srv.get('root_dir') or srv['notebook_dir'])
        return root_dir / path
    raise FileNotFoundError(FILE_ERROR.format('path'))
    