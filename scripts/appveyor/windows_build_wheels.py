from subprocess import check_call
import os
import sys
import shutil

SCRIPT_DIR = os.path.dirname(__file__)
REPO_DIR = os.path.abspath(os.getcwd())
ROOT_DIR = os.path.abspath(os.path.join(REPO_DIR, '..'))

print('ROOT_DIR: %s' % ROOT_DIR)
print('REPO_DIR: %s' % REPO_DIR)

from wheel_build_utils import push_dir, push_env
from windows_build_common import DEFAULT_PY_ENVS, venv_paths

def build_wheels(py_envs=DEFAULT_PY_ENVS):

    # Install Chemfiles
    chfl_build_dir = os.path.join(ROOT_DIR, 'chfl-build')
    os.mkdir(chfl_build_dir)
    chfl_install_dir = os.path.join(ROOT_DIR, 'chfl')
    check_call([
        'cmake', '-DCMAKE_INSTALL_PREFIX:PATH=%s' % chfl_install_dir,
        '-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON',
        '-DCMAKE_BUILD_TYPE:STRING=Release',
        '-G', 'Visual Studio 15 2017 Win64',
        '../chemfiles/'], cwd=chfl_build_dir)
    check_call(['cmake',  '--build',  '.',  '--target', 'install', '--config', 'Release'], cwd=chfl_build_dir)

    for py_env in py_envs:
        python_executable, \
                python_include_dir, \
                python_library, \
                pip, \
                path = venv_paths(py_env)

        with push_env(PATH='%s%s%s' % (path, os.pathsep, os.environ['PATH'])):

            # Install dependencies
            requirements_file = os.path.join(REPO_DIR, 'requirements-dev.txt')
            if os.path.exists(requirements_file):
                check_call([pip, 'install', '--upgrade', '-r', requirements_file])
            check_call([pip, 'install', 'cmake'])
            check_call([pip, 'install', 'scikit_build'])

            build_type = 'Release'

            # Generate wheel
            check_call([
                python_executable,
                'setup.py', 'bdist_wheel',
                '--build-type', build_type, '-G', 'Visual Studio 15 2017 Win64',
                '--',
                '-DPYTHON_EXECUTABLE:FILEPATH=%s' % python_executable,
                '-DPYTHON_INCLUDE_DIR:PATH=%s' % python_include_dir,
                '-DPYTHON_LIBRARY:FILEPATH=%s' % python_library,
                '-DLEMON_EXTERNAL_CHEMFILES:BOOL=ON',
                '-DCHEMFILES_ROOT_DIR:PATH=%s' % os.path.join(chfl_install_dir),
                '-DLEMON_BUILD_PROGS:BOOL=OFF',
            ])
            # Cleanup
            check_call([python_executable, 'setup.py', 'clean'])

if __name__ == '__main__':
    build_wheels()
