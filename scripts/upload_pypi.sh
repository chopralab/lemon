#!/bin/bash
pip install twine
SCRIPT_DIR=$(cd $(dirname $0) || exit 1; pwd)
echo -e "[pypi]" >> ~/.pypirc
echo -e "username = frodofine" >> ~/.pypirc
echo -e "password = $PYPI_PASSWORD" >> ~/.pypirc
twine upload --repository-url https://test.pypi.org/ dist/*